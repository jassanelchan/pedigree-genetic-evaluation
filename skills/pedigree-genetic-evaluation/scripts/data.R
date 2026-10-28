ge_descriptive <- function(data, traits, stage) {
  rows <- lapply(seq_len(nrow(traits)), function(j) {
    value <- data[[paste0("y", j)]]
    value <- value[is.finite(value)]
    n <- length(value)
    average <- if (n) mean(value) else NA_real_
    deviation <- if (n > 1L) sd(value) else NA_real_
    data.frame(trait = traits$name[j], label = traits$label[j], unit = traits$unit[j],
               stage = stage, n = n, mean = average, SD = deviation,
               median = if (n) median(value) else NA_real_,
               min = if (n) min(value) else NA_real_, max = if (n) max(value) else NA_real_,
               CV_percent = if (is.finite(average) && average > 0) deviation / average * 100 else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

ge_prepare_phenotypes <- function(analysis, base) {
  source <- analysis$phenotypes
  raw <- ge_read_tabular(source, base)
  original_rows <- seq_len(nrow(raw))
  original_id <- ge_id(ge_column(raw, source$columns$id, "phenotype animal ID"))
  keep_id <- !is.na(original_id)
  omitted <- sum(!keep_id)
  input_rows <- nrow(raw)
  raw <- raw[keep_id, , drop = FALSE]
  original_rows <- original_rows[keep_id]
  original_id <- original_id[keep_id]
  family <- analysis$model
  if (is.null(family) || !family %in% c("animal", "litter", "repeatability", "general")) {
    ge_stop("Model must be animal, litter, repeatability or general", "unsupported_model")
  }
  if (family != "general" && length(analysis$random)) {
    ge_stop("Named random terms require model general", "unsupported_model")
  }
  if (family != "litter" && length(analysis$shared)) {
    ge_stop("shared is only valid for the litter shortcut; use named random terms in general",
            "inapplicable_shared_definition")
  }
  definitions <- analysis$traits
  if (!length(definitions)) ge_stop("At least one trait is required", "missing_trait_definition")
  traits <- do.call(rbind, lapply(definitions, function(trait) {
    if (is.null(trait$name) || is.null(trait$column) || is.null(trait$unit) ||
        !nzchar(trait$name) || !nzchar(trait$unit)) {
      ge_stop("Each trait needs name, column and verified unit", "missing_trait_definition")
    }
    record_type <- trait$record_type
    if (is.null(record_type)) record_type <- if (family == "repeatability") "repeated" else "single"
    if (!record_type %in% c("single", "repeated")) {
      ge_stop("record_type must be single or repeated", "invalid_record_type")
    }
    data.frame(name = trait$name, label = if (is.null(trait$label)) trait$name else trait$label,
               unit = trait$unit, record_type = record_type, stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(traits$name)) ge_stop("Trait names must be unique")
  nt <- nrow(traits)
  if (family == "general" && length(unique(traits$record_type)) > 1L &&
      any(vapply(definitions, function(trait) is.null(trait$record_filter), logical(1)))) {
    ge_stop("Mixed single and repeated traits need an explicit record_filter for every trait",
            "missing_record_filter")
  }
  data <- data.frame(record_id = original_rows, original_id = original_id,
                     stringsAsFactors = FALSE)
  for (j in seq_len(nt)) {
    definition <- definitions[[j]]
    distribution <- definition$distribution
    if (!is.null(distribution) && !distribution %in% c("gaussian", "count_gaussian")) {
      ge_stop("Only Gaussian and count_gaussian traits are supported", "unsupported_model")
    }
    value <- ge_numeric(ge_column(raw, definition$column, definition$name),
                         definition$name, source$missing_values)
    member <- rep(TRUE, nrow(raw))
    filter <- definition$record_filter
    if (!is.null(filter)) {
      allowed <- ge_id(unlist(filter$values))
      if (!length(allowed) || anyNA(allowed) || anyDuplicated(allowed)) {
        ge_stop("record_filter needs distinct nonmissing values", "invalid_record_filter")
      }
      member <- ge_id(ge_column(raw, filter$column, "trait record filter")) %in% allowed
      if (any(!member & !is.na(value))) {
        ge_stop(paste("Observed trait value outside its declared record_filter:", definition$name),
                "observation_outside_record_filter")
      }
    }
    if (identical(distribution, "count_gaussian") &&
        any(value < 0 | abs(value - round(value)) > 1e-8, na.rm = TRUE)) {
      ge_stop(paste("Count trait must contain nonnegative integers:", definition$name),
              "invalid_count_trait")
    }
    data[[paste0("y", j)]] <- value
    data[[paste0("member", j)]] <- member
  }
  supplement <- NULL
  if (!is.null(source$columns$sire) || !is.null(source$columns$dam)) {
    supplement <- data.frame(
      original_id = original_id,
      sire_original = if (is.null(source$columns$sire)) NA_character_ else
        ge_parent_id(ge_column(raw, source$columns$sire, "supplement sire")),
      dam_original = if (is.null(source$columns$dam)) NA_character_ else
        ge_parent_id(ge_column(raw, source$columns$dam, "supplement dam")),
      stringsAsFactors = FALSE)
  }
  raw_trait_data <- data
  y_columns <- paste0("y", seq_len(nt))
  usable <- rowSums(!is.na(data[y_columns])) > 0L
  missing_all <- sum(!usable)
  data <- data[usable, , drop = FALSE]
  raw <- raw[usable, , drop = FALSE]
  if (!nrow(data)) ge_stop("No usable phenotype records remain", "empty_analysis")
  if (!is.null(analysis$event_date)) {
    data$event_date <- ge_date(ge_column(raw, analysis$event_date, "event date"), source, base)
    data$event_year <- as.integer(format(data$event_date, "%Y"))
  }
  if (!is.null(analysis$birth_date)) {
    data$birth_date <- ge_date(ge_column(raw, analysis$birth_date, "animal birth date"), source, base)
    data$birth_year <- as.integer(format(data$birth_date, "%Y"))
  }
  keys <- analysis$record_key
  if (is.null(keys)) {
    if (any(traits$record_type == "repeated")) {
      ge_stop("Repeated records need an explicit record_key (animal plus event identifier/date)",
              "missing_record_key")
    }
    key_values <- data.frame(id = data$original_id)
  } else {
    if (any(traits$record_type == "repeated") &&
        (length(keys) < 2L || !as.character(source$columns$id) %in% as.character(unlist(keys)))) {
      ge_stop("Repeated record_key must include the mapped animal column and an event",
              "missing_record_key")
    }
    key_values <- as.data.frame(lapply(keys, function(column) {
      if (!is.null(analysis$event_date) && identical(column, analysis$event_date)) return(as.character(data$event_date))
      if (!is.null(analysis$birth_date) && identical(column, analysis$birth_date)) return(as.character(data$birth_date))
      ge_id(ge_column(raw, column, "record key"))
    }))
    if (anyNA(key_values)) ge_stop("Record keys contain missing values", "missing_record_key")
  }
  if (anyDuplicated(key_values)) {
    ge_stop("Duplicate phenotype record keys require resolution", "duplicate_records")
  }
  for (j in seq_len(nt)) if (traits$record_type[j] == "single") {
    observed <- !is.na(data[[y_columns[j]]])
    if (anyDuplicated(data$original_id[observed])) {
      ge_stop(paste("A single-record trait occurs more than once for an animal:", traits$name[j]),
              "duplicate_single_trait")
    }
  }
  before <- data
  qc <- analysis$qc
  method <- if (is.null(qc$method)) "sd" else qc$method
  if (!method %in% c("sd", "none")) {
    ge_stop("QC supports single-pass sd or explicit none", "unsupported_qc")
  }
  multiplier <- if (is.null(qc$sd)) 3 else as.numeric(qc$sd)
  if (length(multiplier) != 1L || !is.finite(multiplier) || multiplier <= 0) {
    ge_stop("QC SD multiplier must be positive")
  }
  group <- if (is.null(qc$group_by)) rep("analysis", nrow(data)) else {
    values <- lapply(qc$group_by, function(column) ge_id(ge_column(raw, column, "QC group")))
    if (any(vapply(values, anyNA, logical(1)))) ge_stop("QC grouping has missing values")
    do.call(paste, c(values, list(sep = "::")))
  }
  audit <- list()
  for (j in seq_len(nt)) for (level in unique(group)) {
    rows <- which(group == level)
    values <- before[[y_columns[j]]][rows]
    finite <- is.finite(values)
    average <- if (any(finite)) mean(values[finite]) else NA_real_
    deviation <- if (sum(finite) > 1L) sd(values[finite]) else NA_real_
    lower <- if (method == "sd") average - multiplier * deviation else -Inf
    upper <- if (method == "sd") average + multiplier * deviation else Inf
    rejected <- if (method == "sd" && is.finite(deviation)) {
      finite & (values < lower | values > upper)
    } else rep(FALSE, length(rows))
    data[[y_columns[j]]][rows[rejected]] <- NA_real_
    audit[[length(audit) + 1L]] <- data.frame(
      trait = traits$name[j], group = level, method = method, n_before = sum(finite),
      mean = average, SD = deviation, lower = lower, upper = upper,
      removed = sum(rejected), stringsAsFactors = FALSE)
  }
  valid_any <- rowSums(!is.na(data[y_columns])) > 0L
  data <- data[valid_any, , drop = FALSE]
  raw <- raw[valid_any, , drop = FALSE]
  if (!nrow(data)) ge_stop("No usable phenotype records remain", "empty_analysis")
  for (j in seq_len(nt)) {
    values <- data[[y_columns[j]]]
    if (sum(!is.na(values)) < 2L || !is.finite(var(values, na.rm = TRUE)) ||
        var(values, na.rm = TRUE) <= 0) {
      ge_stop(paste("Trait has insufficient nonconstant observations:", traits$name[j]),
              "uninformative_trait")
    }
  }
  fixed <- list(list(name = "mean", column = "fx1", type = "class",
                     traits = seq_len(nt), labels = "mean", levels = 1L))
  data$fx1 <- 1L
  definitions_fixed <- analysis$fixed
  fixed_names <- vapply(definitions_fixed, function(effect) {
    if (is.null(effect$name)) "" else effect$name
  }, character(1))
  if (any(!nzchar(fixed_names)) || anyDuplicated(c("mean", fixed_names))) {
    ge_stop("Fixed effects need distinct nonempty names", "invalid_fixed_effect")
  }
  for (index in seq_along(definitions_fixed)) {
    definition <- definitions_fixed[[index]]
    active <- ge_trait_indices(definition$traits, traits)
    needed <- rowSums(!is.na(data[y_columns[active]])) > 0L
    column_name <- paste0("fx", index + 1L)
    type <- if (is.null(definition$type)) "class" else definition$type
    if (!is.null(definition$from_date)) {
      date_name <- if (definition$from_date == "event") "event_date" else "birth_date"
      dates <- data[[date_name]]
      if (is.null(dates)) ge_stop(paste("Date source missing for fixed effect", definition$name))
      value <- switch(definition$grouping,
                      year = format(dates, "%Y"),
                      month = format(dates, "%Y-%m"),
                      quarter = paste(format(dates, "%Y"), paste0("Q", (as.integer(format(dates, "%m")) - 1L) %/% 3L + 1L), sep = "-"),
                      ge_stop("Date grouping must be year, month or quarter"))
      value[is.na(dates)] <- NA_character_
    } else value <- ge_column(raw, definition$column, definition$name)
    if (type == "class") {
      value <- ge_id(value)
      if (anyNA(value[needed])) {
        ge_stop(paste("Fixed effect contains missing categories:", definition$name),
                "missing_fixed_effect")
      }
      levels <- sort(unique(value[needed]), method = "radix")
      codes <- rep(1L, nrow(data))
      codes[needed] <- match(value[needed], levels)
      data[[column_name]] <- codes
      fixed[[length(fixed) + 1L]] <- list(name = definition$name, column = column_name,
                                          type = type, traits = active, labels = levels,
                                          levels = length(levels))
    } else if (type == "cov") {
      values <- rep(0, nrow(data))
      values[needed] <- ge_numeric(value[needed], definition$name, source$missing_values)
      if (anyNA(values[needed])) {
        ge_stop(paste("Covariate contains missing values:", definition$name), "missing_fixed_effect")
      }
      data[[column_name]] <- values
      fixed[[length(fixed) + 1L]] <- list(name = definition$name, column = column_name,
                                          type = type, traits = active)
    } else ge_stop("Fixed effect type must be class or cov")
  }
  random_definitions <- if (family == "general") analysis$random else list()
  shared_role <- NULL
  if (family == "litter") {
    shared_role <- "litter"
    random_definitions <- list(c(list(name = "shared", role = "litter"), analysis$shared))
  }
  if (family == "repeatability") {
    shared_role <- "pe"
    random_definitions <- list(list(name = "shared", role = "pe"))
  }
  random <- list()
  signatures <- character()
  for (index in seq_along(random_definitions)) {
    definition <- random_definitions[[index]]
    name <- definition$name
    role <- definition$role
    if (is.null(name) || length(name) != 1L || !grepl("^[A-Za-z][A-Za-z0-9_]*$", name) ||
        name %in% c("additive", "residual") ||
        (family == "general" && name == "shared") || name %in% names(random)) {
      ge_stop("Random terms need unique safe names, excluding reserved component names",
              "invalid_random_effect")
    }
    if (is.null(role) || !role %in% c("litter", "pe", "iid")) {
      ge_stop("Independent random role must be litter, pe or iid", "unsupported_random_effect")
    }
    active <- ge_trait_indices(definition$traits, traits)
    needed <- rowSums(!is.na(data[y_columns[active]])) > 0L
    column <- if (family == "general") paste0("re", index) else "shared"
    if (role == "pe") {
      if (!is.null(definition$column) || !is.null(definition$dam)) {
        ge_stop("Permanent environment is grouped by animal ID, not a supplied alternative",
                "invalid_pe_group")
      }
      for (trait in active) {
        if (!anyDuplicated(data$original_id[!is.na(data[[y_columns[trait]]])])) {
          ge_stop(paste("No repeated observations identify permanent environment for", traits$name[trait]),
                  "unidentified_pe")
        }
      }
      labels <- data$original_id
    } else if (!is.null(definition$column)) {
      labels <- ge_id(ge_column(raw, definition$column, paste("grouping for", name)))
    } else if (role == "litter" && !is.null(definition$dam) && !is.null(analysis$birth_date)) {
      dam <- ge_parent_id(ge_column(raw, definition$dam, "dam for litter"))
      labels <- paste(dam, data$birth_date, sep = "::")
      labels[is.na(dam) | is.na(data$birth_date)] <- NA_character_
    } else {
      ge_stop(paste("A grouping column or a verified litter definition is required for", name),
              if (role == "litter") "missing_litter_definition" else "missing_random_group")
    }
    if (anyNA(labels[needed])) {
      ge_stop(paste("Missing grouping values for", name),
              if (role == "litter") "missing_litter" else "missing_random_group")
    }
    levels <- sort(unique(labels[needed]), method = "radix")
    codes <- rep(1L, nrow(data))
    codes[needed] <- match(labels[needed], levels)
    data[[column]] <- codes
    partition <- match(labels[needed], unique(labels[needed]))
    signature <- paste(paste(active, collapse = ","), paste(partition, collapse = ","), sep = ":")
    if (signature %in% signatures) {
      ge_stop("Two independent random terms have identical grouping and active traits",
              "confounded_random_effects")
    }
    signatures <- c(signatures, signature)
    random[[name]] <- list(name = name, role = role, relationship = "identity",
                           column = column, traits = active, levels = length(levels),
                           labels = levels)
  }
  for (name in names(random)) {
    term <- random[[name]]
    for (trait in term$traits) {
      observed <- !is.na(data[[y_columns[trait]]])
      groups <- data[[term$column]][observed]
      frequencies <- table(groups)
      if (length(frequencies) < 2L || max(frequencies) < 2L) {
        ge_stop(paste("Random term", name, "has no identifiable replicated groups for",
                      traits$name[trait]), "unidentified_random_effect")
      }
    }
  }
  if (length(random) > 1L) {
    pairs <- combn(names(random), 2L, simplify = FALSE)
    for (pair in pairs) {
      first <- random[[pair[1L]]]
      second <- random[[pair[2L]]]
      for (trait in intersect(first$traits, second$traits)) {
        observed <- !is.na(data[[y_columns[trait]]])
        a <- data[[first$column]][observed]
        b <- data[[second$column]][observed]
        if (identical(match(a, unique(a)), match(b, unique(b)))) {
          ge_stop(paste("Random terms", first$name, "and", second$name,
                        "have identical grouping for", traits$name[trait]),
                  "confounded_random_effects")
        }
      }
    }
  }
  list(data = data, raw_data = raw_trait_data, traits = traits, fixed = fixed,
       random = random, shared_role = shared_role, family = family, supplement = supplement,
       descriptive = rbind(ge_descriptive(before, traits, "before"),
                           ge_descriptive(data, traits, "after")),
       qc_audit = do.call(rbind, audit),
       input_audit = data.frame(
         metric = c("input_rows", "missing_id_rows", "all_traits_missing_rows",
                    "retained_records", "retained_animals"),
         value = c(input_rows, omitted, missing_all, nrow(data),
                   length(unique(data$original_id)))))
}


ge_prepare_pedigree <- function(source, phenotype, targets, base, allow_observed_founders = FALSE) {
  read_parents <- function(table_source) {
    raw <- ge_read_tabular(table_source, base)
    data.frame(original_id = ge_id(ge_column(raw, table_source$columns$id, "pedigree animal ID")),
               sire_original = ge_parent_id(ge_column(raw, table_source$columns$sire, "sire ID")),
               dam_original = ge_parent_id(ge_column(raw, table_source$columns$dam, "dam ID")),
               stringsAsFactors = FALSE)
  }
  ped <- read_parents(source)
  for (extra in source$supplements) ped <- rbind(ped, read_parents(extra))
  ped <- ped[!is.na(ped$original_id), , drop = FALSE]
  if (!is.null(phenotype$supplement)) ped <- rbind(ped, phenotype$supplement)
  merge_parent <- function(value, id) {
    known <- unique(value[!is.na(value)])
    if (length(known) > 1L) ge_stop(paste("Conflicting parents for", id), "conflicting_pedigree")
    if (length(known)) known[[1L]] else NA_character_
  }
  if (anyDuplicated(ped$original_id)) {
    repeated_ids <- unique(ped$original_id[duplicated(ped$original_id)])
    repeated_rows <- which(ped$original_id %in% repeated_ids)
    groups <- split(repeated_rows, ped$original_id[repeated_rows])
    merged <- do.call(rbind, lapply(groups, function(rows) {
      id <- ped$original_id[rows[1L]]
      data.frame(original_id = id,
                 sire_original = merge_parent(ped$sire_original[rows], id),
                 dam_original = merge_parent(ped$dam_original[rows], id), stringsAsFactors = FALSE)
    }))
    ped <- rbind(ped[!ped$original_id %in% repeated_ids, , drop=FALSE], merged)
    rownames(ped) <- NULL
  }
  observed <- unique(phenotype$data$original_id)
  referenced <- unique(c(ped$sire_original, ped$dam_original))
  represented <- unique(c(ped$original_id, referenced[!is.na(referenced)]))
  missing_targets <- setdiff(targets, represented)
  if (length(missing_targets) && !isTRUE(allow_observed_founders)) {
    ge_stop(paste("Prediction IDs absent from supplied pedigree:", paste(head(missing_targets, 5L), collapse = ", ")),
            "prediction_missing_pedigree")
  }
  missing_observed <- setdiff(observed, represented)
  if (length(missing_observed) && !isTRUE(allow_observed_founders)) {
    ge_stop(paste("Observed animals missing from pedigree and without supplied parent records:",
                  paste(head(missing_observed, 5L), collapse = ", ")), "observed_missing_pedigree")
  }
  parent_only <- setdiff(unique(c(ped$sire_original, ped$dam_original, missing_observed, targets)), ped$original_id)
  parent_only <- parent_only[!is.na(parent_only)]
  if (length(parent_only)) {
    ped <- rbind(ped, data.frame(original_id = parent_only, sire_original = NA_character_,
                                 dam_original = NA_character_, stringsAsFactors = FALSE))
  }
  starts <- unique(c(observed, targets))
  selected <- rep(FALSE, nrow(ped))
  frontier <- match(starts, ped$original_id)
  if (anyNA(frontier)) ge_stop("Requested animals not represented in pedigree")
  while (length(frontier)) {
    frontier <- unique(frontier[!selected[frontier]])
    if (!length(frontier)) break
    selected[frontier] <- TRUE
    parents <- c(ped$sire_original[frontier], ped$dam_original[frontier])
    frontier <- match(parents[!is.na(parents)], ped$original_id)
  }
  ped <- ped[selected, , drop = FALSE]
  ped <- ped[order(ped$original_id, method = "radix"), , drop = FALSE]
  rownames(ped) <- NULL
  ped$id <- seq_len(nrow(ped))
  ped$sire <- match(ped$sire_original, ped$original_id)
  ped$dam <- match(ped$dam_original, ped$original_id)
  ped$sire[is.na(ped$sire)] <- 0L
  ped$dam[is.na(ped$dam)] <- 0L
  if (any(ped$id == ped$sire | ped$id == ped$dam)) ge_stop("Self-parent relationship", "invalid_pedigree")
  if (any(ped$sire > 0L & ped$sire == ped$dam)) ge_stop("An animal cannot have the same sire and dam ID", "invalid_pedigree")
  edges_parent <- c(ped$sire[ped$sire > 0L], ped$dam[ped$dam > 0L])
  edges_child <- c(ped$id[ped$sire > 0L], ped$id[ped$dam > 0L])
  children <- split(edges_child, edges_parent)
  indegree <- as.integer(ped$sire > 0L) + as.integer(ped$dam > 0L)
  queue <- integer(nrow(ped))
  founders <- which(indegree == 0L)
  queue[seq_along(founders)] <- founders
  front <- 1L
  back <- length(founders)
  ranking <- integer(nrow(ped))
  visited <- 0L
  while (front <= back) {
    current <- queue[front]
    front <- front + 1L
    visited <- visited + 1L
    ranking[current] <- visited
    offspring <- children[[as.character(current)]]
    if (length(offspring)) {
      indegree[offspring] <- indegree[offspring] - 1L
      ready <- offspring[indegree[offspring] == 0L]
      if (length(ready)) {
        queue[back + seq_along(ready)] <- ready
        back <- back + length(ready)
      }
    }
  }
  if (visited != nrow(ped)) ge_stop("Cycle in the selected pedigree", "pedigree_cycle")
  ped$sort <- ranking
  attr(ped, "parent_only_founders") <- intersect(parent_only, ped$original_id)
  ped[c("id", "sire", "dam", "sort", "original_id")]
}

ge_initial_components <- function(data, nt, shared_role = NULL, random = NULL) {
  variance <- vapply(seq_len(nt), function(j) var(data[[paste0("y", j)]], na.rm = TRUE), numeric(1))
  if (is.null(random)) {
    random <- if (is.null(shared_role)) list() else
      list(shared = list(name = "shared", traits = seq_len(nt)))
  }
  counts <- integer(nt)
  for (term in random) counts[ge_effect_traits(term, nt)] <- counts[ge_effect_traits(term, nt)] + 1L
  matrix_for <- function(fraction, active) {
    value <- matrix(0, nt, nt)
    correlation <- matrix(0.1 / max(1L, length(active) - 1L), length(active), length(active))
    diag(correlation) <- 1
    standard <- sqrt(variance[active] * fraction[active])
    value[active, active] <- correlation * outer(standard, standard)
    value
  }
  result <- list(additive = matrix_for(rep(0.3, nt), seq_len(nt)))
  for (term in random) {
    fraction <- ifelse(counts > 0L, 0.1 / pmax(1L, counts), 0)
    result[[term$name]] <- matrix_for(fraction, ge_effect_traits(term, nt))
  }
  result$residual <- matrix_for(ifelse(counts > 0L, 0.6, 0.7), seq_len(nt))
  result
}

ge_prepare_context <- function(spec, analysis, base, software) {
  inbreeding <- ge_boolean(analysis$inbreeding, "inbreeding", TRUE)
  allow_founders <- ge_boolean(analysis$allow_observed_founders, "allow_observed_founders", FALSE)
  phenotype <- ge_prepare_phenotypes(analysis, base)
  targets <- if (is.null(analysis$predict_ids)) character() else as.character(unlist(analysis$predict_ids))
  pedigree_source <- if (is.null(analysis$pedigree)) spec$pedigree else analysis$pedigree
  if (is.null(pedigree_source)) ge_stop("Pedigree source is required")
  ped <- ge_prepare_pedigree(pedigree_source, phenotype, targets, base,
                             allow_founders)
  data <- phenotype$data
  data$id <- match(data$original_id, ped$original_id)
  if (anyNA(data$id)) ge_stop("Internal phenotype/pedigree ID mismatch")
  nt <- nrow(phenotype$traits)
  random <- phenotype$random
  for (name in names(random)) if (random[[name]]$role == "pe") {
    data[[random[[name]]$column]] <- data$id
    random[[name]]$levels <- nrow(ped)
    random[[name]]$labels <- ped$original_id
  }
  initial <- ge_initial_components(data, nt, random = random)
  all_random <- c(list(additive = list(name = "additive", role = "additive",
                                        relationship = "pedigree", column = "id",
                                        traits = seq_len(nt), levels = nrow(ped),
                                        labels = ped$original_id)), random)
  active_components <- c(lapply(all_random, function(term) term$traits), list(residual = seq_len(nt)))
  constraints <- lapply(active_components, function(active) {
    value <- matrix(TRUE, nt, nt)
    value[active, active] <- FALSE
    value
  })
  for (zero in analysis$zero_covariances) {
    role <- zero$component
    i <- suppressWarnings(as.numeric(unlist(zero$trait1)))
    j <- suppressWarnings(as.numeric(unlist(zero$trait2)))
    if (length(role) != 1L || !role %in% names(initial) || length(i) != 1L ||
        length(j) != 1L || !is.finite(i) || !is.finite(j) ||
        i != floor(i) || j != floor(j) || i < 1L || j < 1L ||
        i > nt || j > nt || i == j) {
      ge_stop("Zero covariance constraints need a known component and two distinct trait indices",
              "unsupported_constraint")
    }
    initial[[role]][i, j] <- initial[[role]][j, i] <- 0
    constraints[[role]][i, j] <- constraints[[role]][j, i] <- TRUE
  }
  for (role in names(analysis$initial)) {
    if (!role %in% names(initial)) ge_stop(paste("Unknown initial component", role))
    rows <- analysis$initial[[role]]
    if (is.matrix(rows)) value <- rows else {
      if (length(rows) != nt || any(lengths(rows) != nt)) {
        ge_stop("Initial matrix must be a full global trait-by-trait nested array",
                "invalid_initial_matrix")
      }
      value <- matrix(as.numeric(unlist(rows)), nt, nt, byrow = TRUE)
    }
    if (!identical(dim(value), c(nt, nt)) || any(!is.finite(value)) ||
        max(abs(value - t(value))) > 1e-10) {
      ge_stop("Initial matrix must be finite and symmetric", "invalid_initial_matrix")
    }
    active <- active_components[[role]]
    inactive <- setdiff(seq_len(nt), active)
    if (length(inactive) && any(value[inactive, ] != 0 | t(value[, inactive, drop = FALSE]) != 0)) {
      ge_stop("Initial values outside an effect's active traits must be zero",
              "invalid_initial_matrix")
    }
    if (any(value == 0 & !constraints[[role]])) {
      ge_stop("Free covariance initial values cannot be zero", "invalid_initial_matrix")
    }
    value[constraints[[role]]] <- 0
    initial[[role]] <- value
  }
  for (role in names(initial)) {
    active <- active_components[[role]]
    if (min(eigen(initial[[role]][active, active, drop = FALSE], symmetric = TRUE,
                  only.values = TRUE)$values) <= 0) {
      ge_stop(paste("Initial active covariance block is not positive definite:", role),
              "invalid_initial_matrix")
    }
  }
  if (nt > 1L) for (i in seq_len(nt - 1L)) for (j in (i + 1L):nt) {
    paired <- !is.na(data[[paste0("y", i)]]) & !is.na(data[[paste0("y", j)]])
    if (!any(paired) && !constraints$residual[i, j]) {
      ge_stop(paste("No real record observes traits", i, "and", j,
                    "together; residual covariance requires an explicit justified zero"),
              "unidentified_residual_covariance")
    }
  }
  for (term in random) if (length(term$traits) > 1L) {
    pairs <- combn(term$traits, 2L, simplify = FALSE)
    for (pair in pairs) {
      first_groups <- data[[term$column]][!is.na(data[[paste0("y", pair[1L])]])]
      second_groups <- data[[term$column]][!is.na(data[[paste0("y", pair[2L])]])]
      if (!length(intersect(first_groups, second_groups)) &&
          !constraints[[term$name]][pair[1L], pair[2L]]) {
        ge_stop(paste("Random term", term$name, "has no groups observed for both traits",
                      pair[1L], "and", pair[2L], "; its covariance needs an explicit justified zero"),
                "unidentified_random_covariance")
      }
    }
  }
  rank_checks <- list()
  coverage <- list()
  for (j in seq_len(nt)) {
    valid <- !is.na(data[[paste0("y", j)]])
    columns <- list()
    for (effect in phenotype$fixed) {
      if (!j %in% ge_effect_traits(effect, nt)) next
      value <- data[[effect$column]][valid]
      if (effect$type == "class") {
        for (level in seq_len(effect$levels)) {
          coverage[[length(coverage) + 1L]] <- data.frame(
            trait = phenotype$traits$name[j], effect = effect$name, level = level,
            label = if (is.null(effect$labels)) as.character(level) else effect$labels[level],
            n = sum(value == level), stringsAsFactors = FALSE)
        }
        value <- factor(value)
        if (nlevels(value) > 1L) columns[[length(columns) + 1L]] <-
          model.matrix(~value)[, -1L, drop = FALSE]
      } else columns[[length(columns) + 1L]] <- matrix(value, ncol = 1L)
    }
    X <- matrix(1, sum(valid), 1L)
    if (length(columns)) X <- do.call(cbind, c(list(X), columns))
    rank <- qr(X)$rank
    rank_checks[[j]] <- data.frame(trait = phenotype$traits$name[j], columns = ncol(X), rank = rank)
    if (rank < ncol(X)) {
      ge_stop(paste("Fixed effect design is rank deficient for", phenotype$traits$name[j]),
              "unidentified_fixed_effects")
    }
  }
  data <- data[order(data$id, data$record_id), , drop = FALSE]
  membership <- as.matrix(data[paste0("member", seq_len(nt))])
  colnames(membership) <- phenotype$traits$name
  rownames(membership) <- NULL
  description <- paste(
    "Species:", analysis$species, "Model:", phenotype$family,
    "Traits:", paste(phenotype$traits$name, collapse = ", "),
    "Fixed:", paste(vapply(phenotype$fixed, function(effect) {
      paste0(effect$name, "[", paste(phenotype$traits$name[effect$traits], collapse = ","), "]")
    }, character(1)), collapse = ", "),
    "Random:", paste(vapply(all_random, function(term) {
      paste0(term$name, ":", term$role, "[",
             paste(phenotype$traits$name[term$traits], collapse = ","), "]")
    }, character(1)), collapse = ", "), "+ residual")
  list(name = analysis$name, species = analysis$species, family = phenotype$family,
       nt = nt, traits = phenotype$traits, data = data, raw_data = phenotype$raw_data,
       fixed = phenotype$fixed, random = all_random, trait_rows = membership,
       shared_role = phenotype$shared_role, ped = ped,
       inbreeding = inbreeding,
       initial = initial, constraints = constraints, software = software,
       timeout = if (is.null(spec$timeout)) 7200 else as.numeric(spec$timeout),
       threads = if (is.null(spec$threads)) 2L else as.integer(spec$threads),
       descriptive = phenotype$descriptive, qc_audit = phenotype$qc_audit,
       input_audit = phenotype$input_audit, rank_checks = do.call(rbind, rank_checks),
       fixed_coverage = do.call(rbind, coverage),
       trend_basis = if (is.null(analysis$trend_basis)) "event_year" else analysis$trend_basis,
       model_description = description)
}
