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
  row_number <- seq_len(nrow(raw))
  original_id <- ge_id(ge_column(raw, source$columns$id, "phenotype animal ID"))
  keep_id <- !is.na(original_id)
  omitted <- sum(!keep_id)
  raw <- raw[keep_id, , drop = FALSE]
  row_number <- row_number[keep_id]
  original_id <- original_id[keep_id]
  definitions <- analysis$traits
  if (!length(definitions) || length(definitions) > 2L) {
    ge_stop("Version 1 supports one or two traits per model", "unsupported_model")
  }
  traits <- do.call(rbind, lapply(definitions, function(t) {
    if (is.null(t$name) || is.null(t$column) || is.null(t$unit)) {
      ge_stop("Each trait needs name, column and verified unit", "missing_trait_definition")
    }
    data.frame(name = t$name, label = if (is.null(t$label)) t$name else t$label,
               unit = t$unit, stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(traits$name)) ge_stop("Trait names must be unique")
  data <- data.frame(record_id = row_number, original_id = original_id,
                     stringsAsFactors = FALSE)
  for (j in seq_len(nrow(traits))) {
    definition <- definitions[[j]]
    if (!is.null(definition$distribution) && !definition$distribution %in% c("gaussian", "count_gaussian")) {
      ge_stop("Binary, ordinal and generalized models are outside version 1", "unsupported_model")
    }
    data[[paste0("y", j)]] <- ge_numeric(ge_column(raw, definition$column, definition$name),
                                         definition$name, source$missing_values)
    if (identical(definition$distribution, "count_gaussian")) {
      value <- data[[paste0("y", j)]]
      if (any(value < 0 | abs(value - round(value)) > 1e-8, na.rm = TRUE)) {
        ge_stop(paste("Count trait must contain nonnegative integers:", definition$name), "invalid_count_trait")
      }
    }
  }
  supplement <- NULL
  if (!is.null(source$columns$sire) || !is.null(source$columns$dam)) {
    supplement <- data.frame(original_id = original_id,
                             sire_original = if (is.null(source$columns$sire)) NA_character_ else ge_parent_id(ge_column(raw, source$columns$sire, "supplement sire")),
                             dam_original = if (is.null(source$columns$dam)) NA_character_ else ge_parent_id(ge_column(raw, source$columns$dam, "supplement dam")),
                             stringsAsFactors = FALSE)
  }
  raw_trait_data <- data
  usable <- rowSums(!is.na(data[paste0("y", seq_len(nrow(traits)))])) > 0L
  missing_all <- sum(!usable)
  data <- data[usable, , drop = FALSE]
  raw <- raw[usable, , drop = FALSE]
  original_id <- original_id[usable]
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
    if (analysis$model == "repeatability") {
      ge_stop("Repeated records need an explicit record_key (animal plus event identifier/date)",
              "missing_record_key")
    }
    key_values <- data.frame(id = original_id)
  } else {
    key_values <- as.data.frame(lapply(keys, function(column) {
      if (!is.null(analysis$event_date) && identical(column, analysis$event_date)) return(as.character(data$event_date))
      if (!is.null(analysis$birth_date) && identical(column, analysis$birth_date)) return(as.character(data$birth_date))
      ge_id(ge_column(raw, column, "record key"))
    }))
    if (anyNA(key_values)) ge_stop("Record keys contain missing values", "missing_record_key")
  }
  if (anyDuplicated(key_values)) ge_stop("Duplicate phenotype record keys require resolution",
                                        "duplicate_records")
  fixed <- list()
  data$fx1 <- 1L
  fixed[[1L]] <- list(name = "mean", column = "fx1", type = "class")
  definitions_fixed <- analysis$fixed
  for (index in seq_along(definitions_fixed)) {
    definition <- definitions_fixed[[index]]
    column_name <- paste0("fx", index + 1L)
    type <- if (is.null(definition$type)) "class" else definition$type
    if (!is.null(definition$from_date)) {
      date_name <- if (definition$from_date == "event") "event_date" else "birth_date"
      dates <- data[[date_name]]
      if (is.null(dates)) ge_stop(paste("Date source missing for fixed effect", definition$name))
      grouping <- definition$grouping
      value <- switch(grouping,
                      year = format(dates, "%Y"),
                      month = format(dates, "%Y-%m"),
                      quarter = paste(format(dates, "%Y"), paste0("Q", (as.integer(format(dates, "%m")) - 1L) %/% 3L + 1L), sep = "-"),
                      ge_stop("Date grouping must be year, month or quarter"))
      value[is.na(dates)] <- NA_character_
    } else value <- ge_column(raw, definition$column, definition$name)
    if (type == "class") {
      value <- ge_id(value)
      if (anyNA(value)) ge_stop(paste("Fixed effect contains missing categories:", definition$name),
                               "missing_fixed_effect")
      levels <- sort(unique(value), method = "radix")
      data[[column_name]] <- match(value, levels)
      fixed[[length(fixed) + 1L]] <- list(name = definition$name, column = column_name,
                                          type = type, labels = levels)
    } else if (type == "cov") {
      value <- ge_numeric(value, definition$name, source$missing_values)
      if (anyNA(value)) ge_stop(paste("Covariate contains missing values:", definition$name),
                               "missing_fixed_effect")
      data[[column_name]] <- value
      fixed[[length(fixed) + 1L]] <- list(name = definition$name, column = column_name, type = type)
    } else ge_stop("Fixed effect type must be class or cov")
  }
  family <- analysis$model
  if (!family %in% c("animal", "litter", "repeatability")) ge_stop("Unsupported model", "unsupported_model")
  shared_role <- NULL
  shared_labels <- NULL
  if (family == "repeatability") {
    if (nrow(traits) != 1L) ge_stop("Multitrait repeated-record models are outside version 1", "unsupported_model")
    shared_role <- "pe"
    shared_labels <- original_id
  }
  if (family == "litter") {
    shared_role <- "litter"
    if (!is.null(analysis$shared$column)) {
      shared_labels <- ge_id(ge_column(raw, analysis$shared$column, "litter ID"))
    } else if (!is.null(analysis$shared$dam) && !is.null(analysis$birth_date)) {
      dam <- ge_parent_id(ge_column(raw, analysis$shared$dam, "dam for litter"))
      shared_labels <- paste(dam, data$birth_date, sep = "::")
      shared_labels[is.na(dam) | is.na(data$birth_date)] <- NA_character_
    } else ge_stop("Litter model needs a verified litter ID or dam plus animal birth date", "missing_litter_definition")
    if (anyNA(shared_labels)) ge_stop("Missing litter IDs cannot be pooled into a shared unknown litter", "missing_litter")
  }
  if (!is.null(shared_role)) {
    shared_levels <- sort(unique(shared_labels), method = "radix")
    data$shared <- match(shared_labels, shared_levels)
  }
  before <- data
  qc <- analysis$qc
  method <- if (is.null(qc$method)) "sd" else qc$method
  if (!method %in% c("sd", "none")) ge_stop("Version 1 provides single-pass sd or explicit none QC", "unsupported_qc")
  multiplier <- if (is.null(qc$sd)) 3 else as.numeric(qc$sd)
  if (!is.finite(multiplier) || multiplier <= 0) ge_stop("QC SD multiplier must be positive")
  group_columns <- qc$group_by
  group <- if (is.null(group_columns)) rep("analysis", nrow(data)) else {
    values <- lapply(group_columns, function(column) ge_id(ge_column(raw, column, "QC group")))
    if (any(vapply(values, anyNA, logical(1)))) ge_stop("QC grouping has missing values")
    do.call(paste, c(values, list(sep = "::")))
  }
  audit <- list()
  for (j in seq_len(nrow(traits))) {
    column <- paste0("y", j)
    for (level in unique(group)) {
      rows <- which(group == level)
      values <- before[[column]][rows]
      finite <- is.finite(values)
      average <- if (any(finite)) mean(values[finite]) else NA_real_
      deviation <- if (sum(finite) > 1L) sd(values[finite]) else NA_real_
      lower <- if (method == "sd") average - multiplier * deviation else -Inf
      upper <- if (method == "sd") average + multiplier * deviation else Inf
      rejected <- if (method == "sd" && is.finite(deviation)) finite & (values < lower | values > upper) else rep(FALSE, length(rows))
      data[[column]][rows[rejected]] <- NA_real_
      audit[[length(audit) + 1L]] <- data.frame(trait = traits$name[j], group = level,
                                              method = method, n_before = sum(finite),
                                              mean = average, SD = deviation, lower = lower,
                                              upper = upper, removed = sum(rejected),
                                              stringsAsFactors = FALSE)
    }
  }
  valid_any <- rowSums(!is.na(data[paste0("y", seq_len(nrow(traits)))])) > 0L
  data <- data[valid_any, , drop = FALSE]
  if (!nrow(data)) ge_stop("No usable phenotype records remain", "empty_analysis")
  for (j in seq_len(nrow(traits))) {
    values <- data[[paste0("y", j)]]
    if (sum(!is.na(values)) < 2L || !is.finite(var(values, na.rm = TRUE)) || var(values, na.rm = TRUE) <= 0) {
      ge_stop(paste("Trait has insufficient nonconstant observations:", traits$name[j]), "uninformative_trait")
    }
  }
  if (family == "repeatability" && !anyDuplicated(data$original_id)) {
    ge_stop("No repeated observations remain to identify permanent environment", "unidentified_pe")
  }
  for (k in seq_along(fixed)) {
    if (fixed[[k]]$type != "class") next
    used <- sort(unique(data[[fixed[[k]]$column]]))
    data[[fixed[[k]]$column]] <- match(data[[fixed[[k]]$column]], used)
    if (!is.null(fixed[[k]]$labels)) fixed[[k]]$labels <- fixed[[k]]$labels[used]
  }
  if (!is.null(shared_role)) data$shared <- match(data$shared, sort(unique(data$shared)))
  list(data = data, raw_data = raw_trait_data, traits = traits, fixed = fixed,
       shared_role = shared_role, family = family, supplement = supplement,
       descriptive = rbind(ge_descriptive(before, traits, "before"), ge_descriptive(data, traits, "after")),
       qc_audit = do.call(rbind, audit),
       input_audit = data.frame(metric = c("input_rows", "missing_id_rows", "all_traits_missing_rows", "retained_records", "retained_animals"),
                                 value = c(sum(keep_id) + omitted, omitted, missing_all, nrow(data), length(unique(data$original_id)))))
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

ge_initial_components <- function(data, nt, shared_role) {
  variance <- vapply(seq_len(nt), function(j) var(data[[paste0("y", j)]], na.rm = TRUE), numeric(1))
  correlation <- matrix(0.1, nt, nt)
  diag(correlation) <- 1
  matrix_for <- function(fraction) {
    standard <- sqrt(variance * fraction)
    correlation * outer(standard, standard)
  }
  result <- list(additive = matrix_for(0.3), residual = matrix_for(if (is.null(shared_role)) 0.7 else 0.6))
  if (!is.null(shared_role)) result$shared <- matrix_for(0.1)
  result
}

ge_prepare_context <- function(spec, analysis, base, software) {
  phenotype <- ge_prepare_phenotypes(analysis, base)
  targets <- if (is.null(analysis$predict_ids)) character() else as.character(unlist(analysis$predict_ids))
  pedigree_source <- if (is.null(analysis$pedigree)) spec$pedigree else analysis$pedigree
  if (is.null(pedigree_source)) ge_stop("Pedigree source is required")
  ped <- ge_prepare_pedigree(pedigree_source, phenotype, targets, base,
                             isTRUE(analysis$allow_observed_founders))
  data <- phenotype$data
  data$id <- match(data$original_id, ped$original_id)
  if (anyNA(data$id)) ge_stop("Internal phenotype/pedigree ID mismatch")
  if (phenotype$family == "repeatability") data$shared <- data$id
  nt <- nrow(phenotype$traits)
  initial <- ge_initial_components(data, nt, phenotype$shared_role)
  constraints <- lapply(initial, function(x) matrix(FALSE, nt, nt))
  if (length(analysis$zero_covariances)) {
    for (zero in analysis$zero_covariances) {
      role <- zero$component
      i <- as.integer(zero$trait1)
      j <- as.integer(zero$trait2)
      if (!role %in% names(initial) || i < 1L || j < 1L || i > nt || j > nt || i == j) {
        ge_stop("Only explicitly declared off-diagonal zero constraints are supported")
      }
      initial[[role]][i,j] <- initial[[role]][j,i] <- 0
      constraints[[role]][i,j] <- constraints[[role]][j,i] <- TRUE
    }
  }
  if (!is.null(analysis$initial)) {
    for (role in names(analysis$initial)) {
      if (!role %in% names(initial)) ge_stop(paste("Unknown initial component", role))
      values <- as.numeric(unlist(analysis$initial[[role]]))
      if (length(values) != nt * nt) ge_stop("Initial matrix has wrong dimensions")
      matrix <- matrix(values, nt, nt, byrow = TRUE)
      if (any(!is.finite(matrix)) || max(abs(matrix - t(matrix))) > 1e-10 ||
          min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
        ge_stop("Initial matrix must be finite, symmetric and positive definite")
      }
      if (any(matrix == 0 & !constraints[[role]])) ge_stop("Free covariance initial values cannot be zero")
      matrix[constraints[[role]]] <- 0
      initial[[role]] <- matrix
    }
  }
  if (nt == 2L && !any(complete.cases(data[c("y1", "y2")])) && !constraints$residual[1L,2L]) {
    ge_stop("No records observe both traits: residual covariance is unidentified; declare a justified zero or change the model",
            "unidentified_residual_covariance")
  }
  rank_checks <- list()
  coverage <- list()
  for (j in seq_len(nt)) {
    valid <- !is.na(data[[paste0("y", j)]])
    columns <- list()
    for (effect in phenotype$fixed) {
      value <- data[[effect$column]][valid]
      if (effect$type == "class") {
        all_levels <- sort(unique(data[[effect$column]]))
        for (level in all_levels) {
          coverage[[length(coverage) + 1L]] <- data.frame(trait = phenotype$traits$name[j],
                                                         effect = effect$name, level = level,
                                                         label = if (is.null(effect$labels)) as.character(level) else effect$labels[level],
                                                         n = sum(value == level))
        }
        value <- factor(value)
        if (nlevels(value) > 1L) columns[[length(columns) + 1L]] <- model.matrix(~value)[,-1L,drop = FALSE]
      } else columns[[length(columns) + 1L]] <- matrix(value, ncol = 1L)
    }
    X <- matrix(1, sum(valid), 1L)
    if (length(columns)) X <- do.call(cbind, c(list(X), columns))
    rank <- qr(X)$rank
    rank_checks[[j]] <- data.frame(trait = phenotype$traits$name[j], columns = ncol(X), rank = rank)
    if (rank < ncol(X)) ge_stop(paste("Fixed effect design is rank deficient for", phenotype$traits$name[j]),
                              "unidentified_fixed_effects")
  }
  data <- data[order(data$id, data$record_id), , drop = FALSE]
  fixed <- phenotype$fixed
  for (k in seq_along(fixed)) if (fixed[[k]]$type == "class") fixed[[k]]$levels <- max(data[[fixed[[k]]$column]])
  ctx <- list(name = analysis$name, species = analysis$species, family = phenotype$family,
              nt = nt, traits = phenotype$traits, data = data, raw_data = phenotype$raw_data,
              fixed = fixed, shared_role = phenotype$shared_role, ped = ped,
              inbreeding = if (is.null(analysis$inbreeding)) TRUE else isTRUE(analysis$inbreeding),
              initial = initial, constraints = constraints, software = software,
              timeout = if (is.null(spec$timeout)) 7200 else as.numeric(spec$timeout),
              threads = if (is.null(spec$threads)) 2L else as.integer(spec$threads),
              descriptive = phenotype$descriptive, qc_audit = phenotype$qc_audit,
              input_audit = phenotype$input_audit,
              rank_checks = do.call(rbind, rank_checks),
              fixed_coverage = do.call(rbind, coverage),
              trend_basis = if (is.null(analysis$trend_basis)) "event_year" else analysis$trend_basis,
              model_description = paste("Species:", analysis$species, "Model:", phenotype$family,
                                          "Fixed:", paste(vapply(fixed, function(x) x$name, character(1)), collapse = ", "),
                                          "Random: additive", if (!is.null(phenotype$shared_role)) paste("+", phenotype$shared_role) else "", "+ residual"))
  ctx
}
