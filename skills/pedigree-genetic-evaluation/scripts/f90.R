# BLUPF90+ adapter. Sourcing this file never starts a program.

f90_number <- function(x) {
  suppressWarnings(as.numeric(gsub("[dD]", "E", x)))
}

f90_assert_directory <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (any(file.exists(file.path(directory, c("renum.par", "renf90.par", "solutions", "blupf90.log"))))) {
    ge_stop("F90 attempt directory already contains program files; use a fresh attempt directory.")
  }
  normalizePath(directory, winslash = "/", mustWork = TRUE)
}

f90_roles <- function(ctx) {
  c("additive", if (!is.null(ctx$shared_role)) "shared", "residual")
}

f90_validate_matrices <- function(ctx, values, fitting = FALSE) {
  for (role in f90_roles(ctx)) {
    x <- values[[role]]
    if (!is.matrix(x) || !identical(dim(x), c(as.integer(ctx$nt), as.integer(ctx$nt))) ||
        any(!is.finite(x)) || max(abs(x - t(x))) > 1e-10 * max(1, max(abs(x)))) {
      ge_stop(paste("Invalid F90 covariance matrix:", role))
    }
    constrained <- ctx$constraints[[role]]
    if (is.null(constrained)) constrained <- matrix(FALSE, ctx$nt, ctx$nt)
    if (!identical(dim(constrained), dim(x)) || any(constrained != t(constrained)) ||
        any(is.na(constrained)) || any(x[constrained] != 0)) {
      ge_stop(paste("Invalid fixed-zero constraint matrix:", role))
    }
    if (any(diag(constrained))) {
      ge_stop("F90 v1 supports fixed-zero covariances, not zero variance components. Remove the component explicitly instead.")
    }
    if (fitting && any(x[!constrained] == 0)) {
      ge_stop(paste("Unconstrained initial zero would be fixed by F90; supply a nonzero start for", role))
    }
    if (min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      ge_stop(paste("F90 needs a positive definite covariance matrix:", role))
    }
  }
  invisible(TRUE)
}

f90_version <- function(lines, kind) {
  pattern <- if (kind == "renum") "RENUMF90 version[[:space:]]+([0-9.]+)" else "BLUPF90\\+ ver\\.[[:space:]]+([0-9.]+)"
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) ge_stop(paste("Cannot verify", kind, "version from program output."), "unsupported_version")
  version <- sub(paste0(".*", pattern, ".*"), "\\1", hit[1])
  expected <- if (kind == "renum") "1.169" else "2.76"
  if (version != expected) {
    ge_stop(paste(kind, "version", version, "is not validated; this adapter is pinned to", expected), "unsupported_version")
  }
  version
}

f90_run_checked <- function(executable, directory, parameter, log_name, ctx) {
  if (is.null(executable) || !file.exists(executable)) ge_stop("Configured F90 executable is missing.", "software_missing")
  input <- file.path(directory, paste0(log_name, ".stdin"))
  writeLines(parameter, input, useBytes = TRUE)
  run <- ge_run(executable, directory = directory, stdin = parameter,
                log_name = log_name, timeout = ctx$timeout,
                env = paste0("OMP_NUM_THREADS=", ctx$threads))
  output <- readLines(file.path(directory, log_name), warn = FALSE)
  if (run$status != 0 || !length(output)) {
    ge_stop(paste("F90 program failed; inspect", log_name, "and its stderr."), "engine_failed")
  }
  list(run = run, output = output)
}

f90_patch_card <- function(card, ctx, values, mode, state) {
  lines <- readLines(card, warn = FALSE)
  key <- trimws(lines)
  residual <- which(key == "RANDOM_RESIDUAL VALUES")
  blocks <- which(key == "(CO)VARIANCES")
  if (length(residual) != 1 || length(blocks) != (1L + as.integer(!is.null(ctx$shared_role)))) {
    ge_stop("Unexpected RENUMF90 covariance blocks; cannot safely build the final card.", "format_mismatch")
  }
  lines[residual + seq_len(ctx$nt)] <- ge_matrix_lines(values$residual)
  lines[blocks[1] + seq_len(ctx$nt)] <- ge_matrix_lines(values$additive)
  if (!is.null(ctx$shared_role)) lines[blocks[2] + seq_len(ctx$nt)] <- ge_matrix_lines(values$shared)
  if (!is.null(state$relationship_file)) {
    key <- trimws(lines)
    groups <- which(key == "RANDOM_GROUP")
    animal <- groups[vapply(groups, function(at) as.integer(trimws(lines[at + 1L])), integer(1)) == state$animal_effect]
    if (length(animal) != 1L) ge_stop("Cannot locate the F90 additive relationship block.", "format_mismatch")
    end <- if (any(groups > animal)) min(groups[groups > animal]) - 1L else length(lines)
    type <- which(key == "RANDOM_TYPE" & seq_along(key) > animal & seq_along(key) <= end)
    file <- which(key == "FILE" & seq_along(key) > animal & seq_along(key) <= end)
    if (length(type) != 1L || length(file) != 1L) ge_stop("Malformed F90 additive relationship block.", "format_mismatch")
    lines[type + 1L] <- "user_file"
    lines[file + 1L] <- state$relationship_file
  }
  lines <- lines[!grepl("^\\s*OPTION\\s+(method|solv_method|sol|se_covar_function|samples_se_covar_function|missing)\\b", lines)]
  options <- c("OPTION missing -999", paste("OPTION method", mode))
  if (mode == "BLUP") options <- c(options, "OPTION solv_method FSPAK", "OPTION sol se", "OPTION residual")
  writeLines(c(lines, options, ""), card, useBytes = TRUE)
  invisible(card)
}

f90_prepare <- function(ctx, directory, values, fitting = TRUE) {
  directory <- f90_assert_directory(directory)
  f90_validate_matrices(ctx, values, fitting)
  nt <- as.integer(ctx$nt)
  if (!nt %in% 1:2 || !length(ctx$fixed)) ge_stop("F90 v1 requires one or two traits and explicit fixed effects.")
  if (any(ctx$ped$id <= 0) || anyDuplicated(ctx$ped$id) ||
      any(!ctx$data$id %in% ctx$ped$id)) ge_stop("Invalid public animal IDs passed to F90.")
  trait_columns <- paste0("y", seq_len(nt))
  input_columns <- c(trait_columns, vapply(ctx$fixed, `[[`, "", "column"), "id",
                     if (!is.null(ctx$shared_role)) "shared")
  if (!all(input_columns %in% names(ctx$data)) || anyDuplicated(input_columns)) ge_stop("F90 input columns are missing or repeated.")
  data <- ctx$data[input_columns]
  for (col in trait_columns) {
    if (any(data[[col]] == -999, na.rm = TRUE)) ge_stop("Observed phenotype conflicts with the F90 missing sentinel -999.")
    data[[col]][is.na(data[[col]])] <- -999
  }
  if (anyNA(data) || any(!is.finite(as.matrix(data)))) ge_stop("Fixed effect or ID is missing/nonfinite in the F90 input.")
  if (any(rowSums(data[trait_columns] != -999) == 0)) ge_stop("F90 input contains an all-traits-missing record.")
  old_digits <- getOption("digits")
  options(digits = 17)
  on.exit(options(digits = old_digits), add = TRUE)
  write.table(data, file.path(directory, "phenotypes.txt"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  ped <- ctx$ped[order(ctx$ped$sort), c("id", "sire", "dam")]
  write.table(ped, file.path(directory, "pedigree.txt"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  fixed_positions <- match(vapply(ctx$fixed, `[[`, "", "column"), input_columns)
  lines <- c("DATAFILE", "phenotypes.txt", "TRAITS", paste(seq_len(nt), collapse = " "),
             "FIELDS_PASSED TO OUTPUT", "", "WEIGHT(S)", "", "RESIDUAL_VARIANCE", ge_matrix_lines(values$residual))
  for (i in seq_along(ctx$fixed)) {
    type <- ctx$fixed[[i]]$type
    if (!type %in% c("class", "cov")) ge_stop("Unsupported F90 fixed effect type.")
    line <- paste(c(rep(fixed_positions[i], nt), if (type == "class") c("cross", "alpha") else "cov"), collapse = " ")
    lines <- c(lines, "EFFECT", line)
  }
  animal_effect <- length(ctx$fixed) + 1L
  shared_effect <- if (!is.null(ctx$shared_role)) animal_effect + 1L else NULL
  lines <- c(lines, "EFFECT", paste(c(rep(match("id", input_columns), nt), "cross", "alpha"), collapse = " "),
             "RANDOM", "animal", "FILE", "pedigree.txt", "FILE_POS", "1 2 3 0 0",
             "PED_DEPTH", "0", "INBREEDING", if (isTRUE(ctx$inbreeding)) "pedigree" else "no-inbreeding",
             "(CO)VARIANCES", ge_matrix_lines(values$additive))
  if (!is.null(ctx$shared_role)) {
    # Export pe as an explicit diagonal effect on the public animal code.
    # RENUM then has the same one-block-per-component structure for litter/pe.
    if (ctx$shared_role == "pe" && any(data$shared != data$id)) ge_stop("PE grouping must equal the animal code.")
    lines <- c(lines, "EFFECT", paste(c(rep(match("shared", input_columns), nt), "cross", "alpha"), collapse = " "),
               "RANDOM", "diagonal", "(CO)VARIANCES", ge_matrix_lines(values$shared))
  }
  lines <- c(lines, "OPTION missing -999", "OPTION ped_search complete", "OPTION alpha_size 40", "")
  writeLines(lines, file.path(directory, "renum.par"), useBytes = TRUE)
  run <- f90_run_checked(ctx$software$renumf90, directory, "renum.par", "renum.stdout.log", ctx)
  version <- f90_version(run$output, "renum")
  if (!file.exists(file.path(directory, "renf90.par")) ||
      !any(grepl("END[[:space:]]+job", run$output))) ge_stop("RENUMF90 did not complete all output files.", "engine_failed")
  state <- list(directory = directory, animal_effect = animal_effect, shared_effect = shared_effect,
                residual_effect = animal_effect + as.integer(!is.null(ctx$shared_role)) + 1L,
                renum_version = version, input_columns = input_columns,
                pedigree_file = sprintf("renadd%02d.ped", animal_effect))
  state$mapping <- f90_read_mapping(ctx, directory, state)
  if (isTRUE(ctx$inbreeding)) {
    relationship <- f90_build_precision_inverse(ctx, directory, state)
    state$relationship_file <- relationship$file
    state$relationship_precision <- relationship$precision
  }
  metadata_keys <- intersect(c("animal_effect", "shared_effect", "residual_effect", "renum_version", "input_columns", "pedigree_file", "relationship_file", "relationship_precision"), names(state))
  ge_write_json(state[metadata_keys],
                file.path(directory, "f90_model_map.json"))
  f90_patch_card(file.path(directory, "renf90.par"), ctx, values, if (fitting) "VCE" else "BLUP", state)
  state
}

f90_read_mapping <- function(ctx, directory, state) {
  path <- file.path(directory, state$pedigree_file)
  if (!file.exists(path)) ge_stop("RENUMF90 did not emit the expected pedigree mapping.", "format_mismatch")
  x <- read.table(path, header = FALSE, colClasses = "character", comment.char = "", quote = "")
  if (ncol(x) != 10) ge_stop("Unsupported RENUMF90 pedigree mapping format; expected the verified 10 columns.", "format_mismatch")
  map <- data.frame(internal_id = f90_number(x[[1]]), id = f90_number(x[[10]]))
  if (anyNA(map) || anyDuplicated(map$id) || anyDuplicated(map$internal_id) ||
      !setequal(map$id, ctx$ped$id)) ge_stop("F90 public ID mapping is incomplete or ambiguous.", "format_mismatch")
  map$original_id <- ctx$ped$original_id[match(map$id, ctx$ped$id)]
  if (anyNA(map$original_id)) ge_stop("F90 mapping lost the original animal identifier.")
  map
}

f90_build_precision_inverse <- function(ctx, directory, state) {
  if (!requireNamespace("Matrix", quietly = TRUE)) ge_stop("F90 pedigree precision export requires the R Matrix package.", "missing_dependency")
  inb_file <- file.path(directory, "renf90.inb")
  if (!file.exists(inb_file)) ge_stop("Higher-precision F90 A-inverse needs verified renf90.inb.", "missing_output")
  inb <- read.table(inb_file, header = FALSE, colClasses = "numeric", comment.char = "", quote = "")
  map <- state$mapping
  n <- nrow(map)
  if (ncol(inb) != 3L || anyNA(inb) || !setequal(inb[[1]], map$id) || anyDuplicated(inb[[1]]) ||
      !setequal(map$internal_id, seq_len(n)) || any(inb[[3]] != map$internal_id[match(inb[[1]], map$id)]) ||
      any(!is.finite(inb[[2]])) || any(inb[[2]] < 0 | inb[[2]] > 1)) {
    ge_stop("Cannot verify F90 F/internal/public ID metadata for precision export.", "format_mismatch")
  }
  F <- rep(NA_real_, n)
  F[inb[[3]]] <- inb[[2]]
  ped <- read.table(file.path(directory, state$pedigree_file), header = FALSE, colClasses = "numeric", comment.char = "", quote = "")
  if (ncol(ped) != 10L || nrow(ped) != n || anyNA(ped) || any(ped[[2]] < 0 | ped[[2]] > n) ||
      any(ped[[3]] < 0 | ped[[3]] > n)) ge_stop("Invalid precision-export pedigree.", "format_mismatch")
  # Verify the engine did not substitute different parents or phantom groups.
  for (position in 2:3) {
    public_parent <- rep(0L, n)
    known <- ped[[position]] > 0L
    public_parent[known] <- map$id[match(ped[[position]][known], map$internal_id)]
    expected <- ctx$ped[[if (position == 2L) "sire" else "dam"]][match(ped[[10]], ctx$ped$id)]
    if (any(public_parent != expected)) ge_stop("F90 renumbered parents differ from the canonical pedigree.", "format_mismatch")
  }
  known_parents <- (ped[[2]] > 0L) + (ped[[3]] > 0L)
  F_sire <- ifelse(ped[[2]] > 0L, F[pmax(1L, ped[[2]])], 0)
  F_dam <- ifelse(ped[[3]] > 0L, F[pmax(1L, ped[[3]])], 0)
  denominator <- 4 - known_parents - F_sire - F_dam
  if (any(!is.finite(denominator)) || any(denominator <= 0)) ge_stop("The pedigree Mendelian sampling variance is nonpositive.", "invalid_pedigree")
  weights <- 4 / denominator
  parts <- lapply(seq_len(n), function(row) {
    ids <- as.integer(unlist(ped[row, 1:3]))
    valid <- ids > 0L
    ids <- ids[valid]
    coefficient <- c(1, -.5, -.5)[valid]
    pairs <- expand.grid(a = seq_along(ids), b = seq_along(ids))
    keep <- ids[pairs$a] <= ids[pairs$b]
    data.frame(i = ids[pairs$a[keep]], j = ids[pairs$b[keep]],
                x = weights[row] * coefficient[pairs$a[keep]] * coefficient[pairs$b[keep]])
  })
  entries <- do.call(rbind, parts)
  inverse <- Matrix::sparseMatrix(i = entries$i, j = entries$j, x = entries$x, dims = c(n, n), symmetric = TRUE)
  triples <- Matrix::summary(inverse)
  lines <- vapply(seq_len(nrow(triples)), function(row)
    paste(triples$i[row], triples$j[row], sprintf("%.17g", triples$x[row])), "")
  file <- "pedigree_Ainv_high_precision.txt"
  writeLines(lines, file.path(directory, file), useBytes = TRUE)
  precision <- list(method = "Henderson D-inverse updates; F90 user_file",
                    F_source = "RENUMF90 1.169 renf90.inb, verified six decimal places",
                    F_is_exact = FALSE, exported_significant_digits = 17L,
                    native_Dinv_code_resolution = .001,
                    max_Dinv_change_from_native_code = max(abs(weights - ped[[4]] / 1000)),
                    n_animals = n, n_upper_triangular_entries = nrow(triples),
                    note = "Same canonical pedigree model; removes native integer Dinv rounding. F retains the program's six-decimal export precision.")
  ge_write_json(precision, file.path(directory, "pedigree_Ainv_precision.json"))
  list(file = file, precision = precision)
}

f90_matrix_block <- function(lines, pattern, nt) {
  hits <- grep(pattern, lines)
  if (!length(hits)) ge_stop(paste("Missing F90 result section:", pattern), "format_mismatch")
  at <- tail(hits, 1)
  rows <- lines[at + seq_len(nt)]
  values <- f90_number(unlist(strsplit(trimws(rows), "[[:space:]]+")))
  if (length(values) != nt * nt || anyNA(values) || any(!is.finite(values))) {
    ge_stop("Malformed F90 final covariance matrix.", "format_mismatch")
  }
  matrix(values, nt, nt, byrow = TRUE)
}

f90_expected_meta <- function(ctx, state) {
  out <- list()
  effects <- c(additive = state$animal_effect, if (!is.null(ctx$shared_role)) c(shared = state$shared_effect), residual = state$residual_effect)
  for (role in names(effects)) {
    constraint <- ctx$constraints[[role]]
    if (is.null(constraint)) constraint <- matrix(FALSE, ctx$nt, ctx$nt)
    for (i in seq_len(ctx$nt)) for (j in i:ctx$nt) if (!constraint[i, j]) {
      out[[length(out) + 1L]] <- data.frame(effect = effects[[role]], role = if (role == "shared") ctx$shared_role else role, trait1 = i, trait2 = j)
    }
  }
  do.call(rbind, out)
}

f90_native_se <- function(ctx, lines) {
  read_at <- function(at) {
    if (length(at) != 1L) return(NULL)
    x <- f90_number(unlist(strsplit(trimws(lines[at + seq_len(ctx$nt)]), "[[:space:]]+")))
    if (length(x) != ctx$nt * ctx$nt || anyNA(x) || any(!is.finite(x)) || any(x < 0)) return(NULL)
    matrix(x, ctx$nt, ctx$nt, byrow = TRUE)
  }
  g <- grep("^[[:space:]]*SE for G[[:space:]]*$", lines)
  n <- 1L + as.integer(!is.null(ctx$shared_role))
  g <- tail(g, n)
  result <- list(additive = if (length(g) == n) read_at(g[1]) else NULL)
  if (!is.null(ctx$shared_role)) result$shared <- if (length(g) == n) read_at(g[2]) else NULL
  result$residual <- read_at(tail(grep("^[[:space:]]*SE for R[[:space:]]*$", lines), 1))
  result
}

f90_read_ai <- function(ctx, directory, state, lines) {
  expected <- f90_expected_meta(ctx, state)
  warning_text <- character()
  path <- file.path(directory, "AI_inv.txt")
  covariance <- NULL
  metadata <- expected[c("role", "trait1", "trait2")]
  checked <- tryCatch({
    if (file.exists(path)) {
      x <- read.table(path, header = FALSE, comment.char = "", quote = "")
      if (ncol(x) != 7 || any(!is.finite(as.matrix(x)))) stop("AI_inv.txt does not match the verified seven-column format")
      n <- max(x[[1]], x[[2]])
      if (n != nrow(expected) || nrow(x) != n * (n + 1) / 2 ||
          any(x[[1]] > x[[2]]) || anyDuplicated(paste(x[[1]], x[[2]]))) stop("AI sampling covariance has missing/duplicate parameters")
      diag_rows <- x[x[[1]] == x[[2]], , drop = FALSE]
      diag_rows <- diag_rows[order(diag_rows[[1]]), , drop = FALSE]
      keys <- paste(diag_rows[[3]], diag_rows[[4]], diag_rows[[5]])
      order_map <- match(keys, paste(expected$effect, expected$trait1, expected$trait2))
      if (anyNA(order_map) || anyDuplicated(order_map)) stop("AI parameter role/trait metadata does not match this model")
      expected_dimension <- ifelse(diag_rows[[3]] == state$residual_effect, 0L, ctx$nt)
      if (any(diag_rows[[6]] != expected_dimension)) stop("AI parameter dimension metadata does not match this model")
      metadata <- expected[order_map, c("role", "trait1", "trait2"), drop = FALSE]
      covariance <- matrix(0, n, n)
      covariance[cbind(x[[1]], x[[2]])] <- x[[7]]
      covariance[cbind(x[[2]], x[[1]])] <- x[[7]]
    } else {
      hits <- grep("inverse of AI matrix \\(Sampling Variance\\)", lines)
      if (!length(hits)) stop("No full inverse AI covariance is available")
      at <- tail(hits, 1)
      n <- nrow(expected)
      rows <- lines[at + seq_len(n)]
      values <- f90_number(unlist(strsplit(trimws(rows), "[[:space:]]+")))
      if (length(values) != n * n || anyNA(values)) stop("Printed inverse AI matrix is incomplete")
      covariance <- matrix(values, n, n, byrow = TRUE)
    }
    if (max(abs(covariance - t(covariance))) > 1e-6 * max(1, max(abs(covariance)))) stop("Inverse AI matrix is not symmetric")
    covariance <- (covariance + t(covariance)) / 2
    ev <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) stop("Inverse AI covariance is not positive definite")
    TRUE
  }, error = function(e) { warning_text <<- conditionMessage(e); FALSE })
  if (!checked) covariance <- NULL
  rownames(metadata) <- NULL
  list(param_cov = covariance, param_meta = metadata, warnings = warning_text)
}

f90_fit <- function(ctx, directory, initial = ctx$initial) {
  state <- f90_prepare(ctx, directory, initial, fitting = TRUE)
  run <- f90_run_checked(ctx$software$blupf90, state$directory, "renf90.par", "vce.stdout.log", ctx)
  version <- f90_version(run$output, "blup")
  if (!any(grepl("FINISHED \\(BLUPF90\\)", run$output))) ge_stop("BLUPF90+ did not finish the variance estimation.", "engine_failed")
  native <- file.path(state$directory, "blupf90.log")
  lines <- if (file.exists(native)) readLines(native, warn = FALSE) else run$output
  components <- list(additive = f90_matrix_block(lines, paste0("Genetic variance\\(s\\) for effect[[:space:]]+", state$animal_effect, "[[:space:]]*$"), ctx$nt))
  if (!is.null(ctx$shared_role)) components$shared <- f90_matrix_block(lines, paste0("Genetic variance\\(s\\) for effect[[:space:]]+", state$shared_effect, "[[:space:]]*$"), ctx$nt)
  components$residual <- f90_matrix_block(lines, "^[[:space:]]*Residual variance\\(s\\)", ctx$nt)
  ai <- f90_read_ai(ctx, state$directory, state, lines)
  round_lines <- grep("In round.*convergence", lines, value = TRUE)
  rounds <- if (length(round_lines)) as.integer(sub(".*In round[[:space:]]+([0-9]+).*", "\\1", tail(round_lines, 1))) else NA_integer_
  convergence <- if (length(round_lines)) f90_number(sub(".*convergence[[:space:]]*=[[:space:]]*([^[:space:]]+).*", "\\1", tail(round_lines, 1))) else NA_real_
  ll <- grep("-2logL[[:space:]]*=", lines, value = TRUE)
  minus2logL <- if (length(ll)) f90_number(sub(".*-2logL[[:space:]]*=[[:space:]]*([^[:space:]]+).*", "\\1", tail(ll, 1))) else NA_real_
  warnings <- c(ai$warnings, grep("NaN|bending|not positive|singular|maximum.*round|maximum.*iteration", run$output, value = TRUE, ignore.case = TRUE))
  for (role in names(components)) {
    ev <- eigen(components[[role]], symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) warnings <- c(warnings, paste("Final", role, "covariance matrix is not positive definite; fixed-VC BLUP will not run."))
    else if (min(ev) / max(ev) < 1e-8) warnings <- c(warnings, paste("Final", role, "covariance matrix is close to a boundary."))
  }
  if (is.na(convergence) || convergence > 1e-8) warnings <- c(warnings, "Final variance convergence criterion is missing or above 1e-8.")
  converged <- is.finite(convergence) && convergence <= 1e-8 &&
    !any(grepl("maximum.*round|maximum.*iteration", warnings, ignore.case = TRUE))
  native_se <- f90_native_se(ctx, lines)
  state$blup_version <- version
  ge_write_json(list(engine = "f90", version = list(renumf90 = state$renum_version, blupf90 = version),
                     rounds = rounds, criterion = convergence, converged = converged, minus2logL = minus2logL,
                     param_meta = ai$param_meta, warnings = unique(warnings)), file.path(state$directory, "fit_status.json"))
  list(engine = "f90", components = components, param_meta = ai$param_meta,
       param_cov = ai$param_cov, converged = converged, rounds = rounds,
       minus2logL = minus2logL, native_se = native_se, warnings = unique(warnings),
       files = list.files(state$directory, full.names = TRUE),
       version = list(renumf90 = state$renum_version, blupf90 = version), state = state)
}

f90_read_fitted <- function(ctx, directory, state) {
  file <- file.path(directory, "yhat_residual")
  if (!file.exists(file)) return(list(data = NULL, note = "Native fitted values were not exported."))
  tryCatch({
    native <- read.table(file, header = FALSE, comment.char = "", quote = "")
    encoded <- read.table(file.path(directory, "renf90.dat"), header = FALSE, comment.char = "", quote = "")
    if (ncol(native) != 2L * ctx$nt || nrow(native) != nrow(ctx$data) ||
        nrow(encoded) != nrow(ctx$data) || ncol(encoded) != length(state$input_columns)) {
      stop("Native fitted-value output has an unsupported shape.")
    }
    internal_id <- encoded[[match("id", state$input_columns)]]
    public_id <- state$mapping$id[match(internal_id, state$mapping$internal_id)]
    if (anyNA(public_id) || any(public_id != ctx$data$id)) stop("Native fitted-value record order cannot be verified.")
    result <- list()
    for (trait in seq_len(ctx$nt)) {
      y <- ctx$data[[paste0("y", trait)]]
      observed <- !is.na(y)
      predicted <- f90_number(native[[trait]])
      residual <- f90_number(native[[ctx$nt + trait]])
      if (anyNA(predicted) || anyNA(residual) || any(!is.finite(predicted)) || any(!is.finite(residual)) ||
          any(!observed != (residual == -999))) stop("Native fitted-value missing masks cannot be verified.")
      tolerance <- 5e-6 * max(1, abs(y[observed]))
      if (any(abs(y[observed] - predicted[observed] - residual[observed]) > tolerance) ||
          any(abs(y[observed] - encoded[[trait]][observed]) > 1e-10 * max(1, abs(y[observed])))) {
        stop("Native fitted values do not reproduce observed phenotypes with their residuals.")
      }
      # F90 writes zero yhat at missing observations; that zero is not a prediction.
      predicted[!observed] <- NA_real_
      result[[trait]] <- data.frame(record_id = ctx$data$record_id, trait = trait, predicted = predicted)
    }
    list(data = do.call(rbind, result), note = NULL)
  }, error = function(e) list(data = NULL, note = conditionMessage(e)))
}

f90_blup <- function(ctx, fit, directory, components = fit$components) {
  f90_validate_matrices(ctx, components, fitting = FALSE)
  directory <- f90_assert_directory(directory)
  if (is.null(fit$state) || is.null(fit$state$mapping)) {
    state <- f90_prepare(ctx, directory, components, fitting = FALSE)
  } else {
    state <- fit$state
    if (isTRUE(ctx$inbreeding) && is.null(state$relationship_file)) {
      ge_stop("Cached F90 fit uses the former integer pedigree precision; rerun REML with the current adapter before BLUP.", "stale_fit_precision")
    }
    needed <- c("renf90.par", "renf90.dat", state$pedigree_file, state$relationship_file)
    optional <- c("renf90.inb", "f90_model_map.json", "pedigree_Ainv_precision.json")
    source_dir <- state$directory
    files <- c(needed, optional[file.exists(file.path(source_dir, optional))])
    if (!all(file.copy(file.path(source_dir, files), directory, overwrite = FALSE))) ge_stop("Cannot stage F90 BLUP input files.")
    state$directory <- directory
    state$mapping <- f90_read_mapping(ctx, directory, state)
    f90_patch_card(file.path(directory, "renf90.par"), ctx, components, "BLUP", state)
  }
  run <- f90_run_checked(ctx$software$blupf90, directory, "renf90.par", "blup.stdout.log", ctx)
  version <- f90_version(run$output, "blup")
  if (!any(grepl("FINISHED \\(BLUPF90\\)", run$output)) ||
      !any(grepl("solutions and s.e. stored", run$output, fixed = TRUE))) ge_stop("F90 BLUP did not finish or did not export prediction errors.", "engine_failed")
  path <- file.path(directory, "solutions")
  x <- read.table(path, header = FALSE, skip = 1, comment.char = "", quote = "")
  if (ncol(x) != 5) ge_stop("Unsupported F90 solutions format; expected trait/effect/level/solution/se.", "format_mismatch")
  x <- x[x[[2]] == state$animal_effect, , drop = FALSE]
  mapping <- state$mapping
  at <- match(x[[3]], mapping$internal_id)
  observation_keys <- unlist(lapply(seq_len(ctx$nt), function(trait) {
    paste(unique(ctx$data$id[!is.na(ctx$data[[paste0("y", trait)]])]), trait, sep = ":")
  }), use.names = FALSE)
  ebv <- data.frame(id = mapping$id[at], original_id = mapping$original_id[at],
                    trait = as.integer(x[[1]]), EBV = f90_number(x[[4]]), SEP = f90_number(x[[5]]),
                    F = NA_real_, has_record = paste(mapping$id[at], x[[1]], sep = ":") %in% observation_keys,
                    stringsAsFactors = FALSE)
  expected <- as.vector(outer(ctx$ped$id, seq_len(ctx$nt), paste, sep = ":"))
  found <- paste(ebv$id, ebv$trait, sep = ":")
  if (anyNA(ebv[c("id", "original_id", "trait", "EBV", "SEP")]) || anyDuplicated(found) ||
      !setequal(found, expected) || any(!is.finite(ebv$EBV)) || any(!is.finite(ebv$SEP)) || any(ebv$SEP < 0)) {
    ge_stop("F90 EBV output does not cover every public animal-by-trait pair with valid SEP.", "format_mismatch")
  }
  inb_file <- file.path(directory, "renf90.inb")
  inb_reason <- NULL
  if (isTRUE(ctx$inbreeding)) {
    if (!file.exists(inb_file)) {
      inb_reason <- "RENUMF90 did not export renf90.inb; F remains unavailable."
    } else {
      inb <- read.table(inb_file, header = FALSE, colClasses = "character", comment.char = "", quote = "")
      if (ncol(inb) != 3) ge_stop("Unsupported renf90.inb format.", "format_mismatch")
      pub <- f90_number(inb[[1]])
      internal <- f90_number(inb[[3]])
      match_map <- match(pub, mapping$id)
      F <- f90_number(inb[[2]])
      if (anyNA(match_map) || anyDuplicated(pub) || !setequal(pub, mapping$id) ||
          any(mapping$internal_id[match_map] != internal) || any(!is.finite(F)) || any(F < 0 | F > 1)) {
        ge_stop("Inbreeding output does not match the verified public/internal ID mapping.", "format_mismatch")
      }
      ebv$F <- F[match(ebv$id, pub)]
    }
  } else inb_reason <- "Inbreeding correction was disabled; no verified pedigree F is returned."
  fitted <- f90_read_fitted(ctx, directory, state)
  ge_write_json(list(engine = "f90", version = version, n_animals = nrow(mapping), nt = ctx$nt,
                     coverage = "complete", F_note = inb_reason, fitted_note = fitted$note,
                     SEP_note = "FSPAK OPTION sol se: prediction standard error"),
                file.path(directory, "blup_status.json"))
  ebv <- ebv[order(ebv$id, ebv$trait), , drop = FALSE]
  rownames(ebv) <- NULL
  attr(ebv, "diagnostics") <- list(F_note = inb_reason, fitted_note = fitted$note, version = version)
  attr(ebv, "fitted") <- fitted$data
  ebv
}
