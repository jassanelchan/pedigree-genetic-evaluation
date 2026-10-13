# DMU v6 release 5.6 pedigree backend. Sourcing this module performs no work.

dmu_number <- function(x) suppressWarnings(as.numeric(gsub("[dD]", "E", x)))

dmu_numeric_rows <- function(path, columns) {
  if (!file.exists(path)) ge_stop(paste("DMU output missing:", path), "missing_output")
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  values <- lapply(strsplit(lines, "[[:space:]]+"), dmu_number)
  if (!length(values) || any(lengths(values) != columns) ||
      any(!is.finite(unlist(values)))) {
    ge_stop(paste("Invalid DMU numeric output:", path), "invalid_output")
  }
  do.call(rbind, values)
}

dmu_runtime <- function(ctx) {
  directory <- ctx$software$dmu_dir
  if (is.null(directory) || !is.character(directory) || length(directory) != 1L) {
    ge_stop("DMU requires software.dmu_dir pointing to the verified v6 R5.6 bin directory",
            "missing_engine")
  }
  directory <- normalizePath(directory, winslash = "/", mustWork = FALSE)
  extension <- if (.Platform$OS.type == "windows") ".exe" else ""
  required <- paste0(c("dmu1", "dmuai", "dmu4"), extension)
  if (!all(file.exists(file.path(directory, required))) &&
      all(file.exists(file.path(directory, "bin", required)))) {
    directory <- file.path(directory, "bin")
  }
  if (!all(file.exists(file.path(directory, required)))) {
    ge_stop(paste("DMU v6 R5.6 executables not found in", directory), "missing_engine")
  }
  list(directory = directory, dmu1 = file.path(directory, required[1L]),
       dmuai = file.path(directory, required[2L]), dmu4 = file.path(directory, required[3L]))
}

dmu_roles <- function(ctx) {
  if (is.null(ctx$shared_role)) c("additive", "residual") else
    c("additive", ctx$shared_role, "residual")
}

dmu_block_key <- function(role) if (role %in% c("litter", "pe")) "shared" else role

dmu_validate_context <- function(ctx, matrices) {
  if (!ctx$nt %in% c(1L, 2L) || nrow(ctx$traits) != ctx$nt) {
    ge_stop("DMU v1 supports one or two explicitly ordered traits", "unsupported_model")
  }
  if (!is.null(ctx$shared_role) && !ctx$shared_role %in% c("litter", "pe")) {
    ge_stop("Unsupported DMU shared relationship", "unsupported_model")
  }
  if (identical(ctx$shared_role, "pe") && ctx$nt != 1L) {
    ge_stop("The current DMU repeatability adapter supports one trait", "unsupported_model")
  }
  if (!is.data.frame(ctx$ped) || !nrow(ctx$ped) || anyDuplicated(ctx$ped$id) ||
      anyNA(ctx$ped[c("id", "sire", "dam", "sort")])) {
    ge_stop("DMU needs a complete unique integer pedigree", "invalid_pedigree")
  }
  ints <- as.matrix(ctx$ped[c("id", "sire", "dam", "sort")])
  if (any(!is.finite(ints)) || any(ints != floor(ints)) || any(ints < 0) ||
      any(ints > .Machine$integer.max) || any(ctx$ped$id == 0)) {
    ge_stop("DMU pedigree codes must fit INTEGER*4; UPG codes are unsupported in v1",
            "unsupported_pedigree")
  }
  if (!all(ctx$data$id %in% ctx$ped$id)) ge_stop("DMU phenotype IDs are absent from pedigree")
  for (parent in c("sire", "dam")) {
    index <- match(ctx$ped[[parent]], ctx$ped$id)
    known <- ctx$ped[[parent]] != 0
    if (anyNA(index[known]) || any(ctx$ped$sort[index[known]] >= ctx$ped$sort[known])) {
      ge_stop("DMU pedigree sort must place every known parent before its offspring",
              "invalid_pedigree")
    }
  }
  if (anyDuplicated(ctx$ped$sort) || any(ctx$ped$sort <= 0)) {
    ge_stop("DMU needs a positive unique validated pedigree sort", "invalid_pedigree")
  }
  for (fx in ctx$fixed) {
    if (!fx$type %in% c("class", "cov") || !fx$column %in% names(ctx$data)) {
      ge_stop("Invalid DMU fixed-effect specification", "unsupported_model")
    }
    v <- ctx$data[[fx$column]]
    if (anyNA(v) || any(!is.finite(v))) ge_stop("DMU fixed effects cannot be missing")
    if (fx$type == "class" &&
        any(v <= 0 | v != floor(v) | v > .Machine$integer.max)) {
      ge_stop("DMU class codes must be positive INTEGER*4 values")
    }
  }
  if (!is.null(ctx$shared_role)) {
    shared <- ctx$data$shared
    if (is.null(shared) || anyNA(shared) ||
        any(shared <= 0 | shared != floor(shared) | shared > .Machine$integer.max)) {
      ge_stop("DMU shared-effect codes must be positive INTEGER*4 values")
    }
    if (identical(ctx$shared_role, "pe")) {
      pairs <- unique(ctx$data[c("id", "shared")])
      if (anyDuplicated(pairs$id) || anyDuplicated(pairs$shared)) {
        ge_stop("Permanent environment must identify exactly one effect per animal")
      }
    }
  }
  for (role in dmu_roles(ctx)) {
    key <- dmu_block_key(role)
    x <- matrices[[key]]
    if (is.null(x) || !identical(dim(x), c(as.integer(ctx$nt), as.integer(ctx$nt))) ||
        any(!is.finite(x)) || max(abs(x - t(x))) > 1e-10 * max(1, max(abs(x))) ||
        any(diag(x) <= 0) || min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      ge_stop(paste("DMU covariance matrix must be finite, symmetric and positive definite:", key),
              "invalid_components")
    }
    restriction <- ctx$constraints[[key]]
    if (!is.null(restriction)) {
      if (!identical(dim(restriction), dim(x)) || anyNA(restriction) ||
          !identical(as.logical(restriction), as.logical(t(restriction))) ||
          any(diag(restriction)) || any(x[restriction] != 0)) {
        ge_stop(paste("Invalid fixed-zero covariance specification:", key), "invalid_constraint")
      }
    }
  }
  invisible(TRUE)
}

dmu_is_fixed <- function(ctx, role, i, j) {
  x <- ctx$constraints[[dmu_block_key(role)]]
  !is.null(x) && isTRUE(x[i, j])
}

dmu_stage <- function(ctx, directory, matrices, task) {
  dmu_validate_context(ctx, matrices)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  filenames <- c("analysis.DIR", "data.txt", "pedigree.txt", "MODINF", "PAROUT", "SOL")
  if (any(file.exists(file.path(directory, filenames)))) {
    ge_stop("DMU needs a fresh independent attempt directory; existing outputs are not overwritten",
            "nonempty_attempt_directory")
  }
  class_effects <- Filter(function(x) identical(x$type, "class"), ctx$fixed)
  cov_effects <- Filter(function(x) identical(x$type, "cov"), ctx$fixed)
  integer_data <- data.frame(animal = ctx$data$id)
  if (length(class_effects)) {
    for (i in seq_along(class_effects)) integer_data[[paste0("fx", i)]] <-
      ctx$data[[class_effects[[i]]$column]]
  }
  shared_index <- NULL
  if (!is.null(ctx$shared_role)) {
    integer_data$shared <- ctx$data$shared
    shared_index <- ncol(integer_data)
  }
  real_data <- as.data.frame(ctx$data[paste0("y", seq_len(ctx$nt))])
  if (length(cov_effects)) {
    for (i in seq_along(cov_effects)) real_data[[paste0("cv", i)]] <-
      ctx$data[[cov_effects[[i]]$column]]
  }
  if (any(rowSums(!is.na(real_data[seq_len(ctx$nt)])) == 0)) {
    ge_stop("DMU was given records with no valid trait observation", "invalid_context")
  }
  finite <- unlist(real_data, use.names = FALSE)
  finite <- finite[is.finite(finite)]
  if (!length(finite) || any(abs(finite) > 1e35)) {
    ge_stop("DMU REAL*4 data values are out of the supported range", "unsupported_data_scale")
  }
  low <- min(finite)
  missing_limit <- low - max(1000, abs(low) * 0.01)
  sentinel <- missing_limit - max(1000, abs(missing_limit) * 0.01)
  real_data[is.na(real_data)] <- sentinel
  combined <- cbind(integer_data, real_data)
  records <- apply(combined, 1L, function(x) paste(sprintf("%.17g", x), collapse = " "))
  writeLines(records, file.path(directory, "data.txt"), useBytes = TRUE)
  sorted_ped <- ctx$ped[order(ctx$ped$sort), c("id", "sire", "dam", "sort")]
  write.table(sorted_ped, file.path(directory, "pedigree.txt"), row.names = FALSE,
              col.names = FALSE, quote = FALSE, sep = " ")

  roles <- dmu_roles(ctx)
  random_cols <- c(1L, shared_index)
  model_cols <- c(if (length(class_effects)) seq_len(length(class_effects)) + 1L, random_cols)
  nrandom <- length(random_cols)
  model <- c(as.character(ctx$nt), rep("0", ctx$nt))
  for (trait in seq_len(ctx$nt)) {
    model <- c(model, paste(c(trait, 0L, length(model_cols), model_cols), collapse = " "))
  }
  model <- c(model, rep(paste(c(nrandom, seq_len(nrandom)), collapse = " "), ctx$nt))
  regression <- paste(c(length(cov_effects), if (length(cov_effects))
    seq_len(length(cov_effects)) + ctx$nt), collapse = " ")
  model <- c(model, rep(regression, ctx$nt), "0")
  priors <- character()
  restrictions <- character()
  parameter_keys <- list()
  for (factor in seq_along(roles)) {
    role <- roles[factor]
    matrix <- matrices[[dmu_block_key(role)]]
    for (row in seq_len(ctx$nt)) for (column in seq_len(row)) {
      priors <- c(priors, paste(c(factor, row, column, sprintf("%.17g", matrix[row, column])),
                               collapse = " "))
      parameter_keys[[length(parameter_keys) + 1L]] <-
        data.frame(factor = factor, role = role, trait1 = row, trait2 = column,
                   fixed = dmu_is_fixed(ctx, role, row, column))
      if (dmu_is_fixed(ctx, role, row, column)) {
        restrictions <- c(restrictions, paste("$VAR_REST COV", factor, "E", paste(row, column, sep = ",")))
      }
    }
  }
  method <- if (task == 1L) if (isTRUE(ctx$stable_reml)) 4L else 1L else 9L
  card <- c("$COMMENT", "Pedigree evaluation; reproducible engine adapter", "",
            paste("$ANALYSE", task, method, "0 0"),
            paste0("$DATA ASCII (", ncol(integer_data), ",", ncol(real_data), ",",
                   sprintf("%.17g", missing_limit), ") data.txt"),
            "$VARIABLE", paste(names(integer_data), collapse = " "),
            paste(names(real_data), collapse = " "), "$MODEL", model,
            paste("$VAR_STR 1 PED", if (isTRUE(ctx$inbreeding)) 1L else 2L,
                  "ASCII pedigree.txt"), restrictions, "$PRIOR", priors)
  if (task == 1L) card <- c(card, "$DMUAI", "10", "1.0d-7", "1.0d-6", "0", "0", "0") else
    card <- c(card, "$RESIDUALS ASCII")
  writeLines(card, file.path(directory, "analysis.DIR"), useBytes = TRUE)
  manifest <- list(engine = "DMU", supported_release = "v6 R5.6", task = task, method = method,
                   integer_columns = names(integer_data), real_columns = names(real_data),
                   expected_observations = as.integer(colSums(!is.na(ctx$data[paste0("y", seq_len(ctx$nt))]))),
                   record_count = nrow(ctx$data), pedigree_count = nrow(ctx$ped),
                   missing_limit = missing_limit, missing_sentinel = sentinel,
                   pedigree_method = if (isTRUE(ctx$inbreeding)) 1L else 2L,
                   roles = roles, covariance_keys = do.call(rbind, parameter_keys))
  ge_write_json(manifest, file.path(directory, "input_manifest.json"))
  list(directory = directory, card = card, manifest = manifest)
}

dmu_check_observations <- function(ctx, text) {
  marker <- grep("DEPENDENT VARIABLES \\(TRAITS\\)", text)
  if (!length(marker)) ge_stop("DMU observation summary is missing", "invalid_output")
  start <- marker[length(marker)] + 1L
  end_markers <- grep("COVARIABLES|No of obs for each trait combination|MODEL STATISTICS", text)
  end_markers <- end_markers[end_markers > start]
  end <- if (length(end_markers)) min(end_markers) - 1L else min(length(text), start + ctx$nt + 6L)
  area <- text[seq.int(start, end)]
  rows <- lapply(seq_len(ctx$nt), function(t) {
    line <- grep(paste0("^[[:space:]]*y", t, "[[:space:]]+"), area, value = TRUE)
    if (length(line) != 1L) ge_stop("DMU trait-count summary is ambiguous", "invalid_output")
    fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    dmu_number(fields[3L])
  })
  actual <- as.integer(unlist(rows))
  expected <- as.integer(colSums(!is.na(ctx$data[paste0("y", seq_len(ctx$nt))])))
  if (!identical(actual, expected)) {
    ge_stop(paste("DMU read different effective trait counts; expected",
                  paste(expected, collapse = ","), "but found", paste(actual, collapse = ",")),
            "observation_mask_mismatch")
  }
  actual
}

dmu_run_backend <- function(ctx, stage, mode) {
  runtime <- dmu_runtime(ctx)
  timeout <- if (is.null(ctx$timeout)) 7200 else as.numeric(ctx$timeout)
  if (!is.finite(timeout) || timeout <= 0) ge_stop("A positive DMU process timeout is required")
  started <- Sys.time()
  remaining <- function() {
    value <- timeout - as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (value <= 0) ge_stop("DMU total attempt timeout exhausted", "engine_timeout")
    value
  }
  threads <- if (is.null(ctx$threads)) 2L else as.integer(ctx$threads)
  if (is.na(threads) || threads < 1L) ge_stop("DMU threads must be a positive integer")
  env <- c(paste0("MKL_NUM_THREADS=", threads), paste0("OMP_NUM_THREADS=", threads))
  first <- ge_run(runtime$dmu1, directory = stage$directory, stdin = stage$card,
                  log_name = "dmu1.stdout.log", timeout = remaining(), env = env)
  if (identical(first$status, 124L)) ge_stop("DMU1 exceeded its attempt timeout", "engine_timeout")
  first_text <- if (file.exists(first$stdout)) readLines(first$stdout, warn = FALSE) else character()
  semantic_errors <- any(grepl("^[[:space:]]*ERRORS?[[:space:]:-]", first_text))
  if (first$status != 0L || semantic_errors || !file.exists(file.path(stage$directory, "MODINF")) ||
      !any(grepl("DMU1: Finished with exit code 0", first_text, fixed = TRUE))) {
    ge_stop(paste("DMU1 preparation failed; inspect", first$stdout, "and", first$stderr),
            "engine_preparation_failed")
  }
  version <- grep("Version:.*Release:", first_text, value = TRUE)
  if (!length(version) || !any(grepl("Version:[[:space:]]*6,[[:space:]]*Release:[[:space:]]*5\\.6([^0-9]|$)", version))) {
    ge_stop("The DMU adapter supports verified version 6 release 5.6 only", "unsupported_engine_version")
  }
  dmu_check_observations(ctx, first_text)
  if (as.numeric(difftime(Sys.time(), started, units = "secs")) >= timeout) {
    ge_stop("DMU total attempt timeout exhausted after preprocessing", "engine_timeout")
  }
  program <- if (mode == "fit") runtime$dmuai else runtime$dmu4
  module <- if (mode == "fit") "DMUAI" else "DMU4"
  second <- ge_run(program, directory = stage$directory,
                   log_name = if (mode == "fit") "dmuai.stdout.log" else "dmu4.stdout.log",
                   timeout = remaining(), env = env)
  if (identical(second$status, 124L)) ge_stop(paste(module, "exceeded its attempt timeout"), "engine_timeout")
  second_text <- if (file.exists(second$stdout)) readLines(second$stdout, warn = FALSE) else character()
  ended <- if (mode == "fit") any(grepl("DMUAI: Finished with exit code 0", second_text, fixed = TRUE)) else
    any(grepl("^[[:space:]]*End of DMU4[[:space:]]*$", second_text))
  if (second$status != 0L || !ended) {
    ge_stop(paste(module, "failed; inspect", second$stdout, "and", second$stderr), "engine_run_failed")
  }
  list(runtime = runtime, runs = list(first, second), preparation = first_text,
       output = second_text, version = trimws(version[1L]), elapsed =
         as.numeric(difftime(Sys.time(), started, units = "secs")))
}

dmu_parse_components <- function(ctx, path) {
  raw <- dmu_numeric_rows(path, 4L)
  if (any(raw[, 1:3] != floor(raw[, 1:3]))) ge_stop("Invalid DMU component indices", "invalid_output")
  roles <- dmu_roles(ctx)
  if (any(raw[, 1L] < 1L | raw[, 1L] > length(roles)) ||
      any(raw[, 2:3] < 1L | raw[, 2:3] > ctx$nt)) {
    ge_stop("DMU covariance output does not match the configured model", "invalid_output")
  }
  row <- pmax(raw[, 2L], raw[, 3L])
  column <- pmin(raw[, 2L], raw[, 3L])
  keys <- paste(raw[, 1L], row, column, sep = ":")
  if (anyDuplicated(keys)) ge_stop("Duplicated DMU covariance entries", "invalid_output")
  components <- list()
  metadata <- data.frame()
  for (factor in seq_along(roles)) {
    role <- roles[factor]
    matrix <- matrix(0, ctx$nt, ctx$nt)
    for (r in seq_len(ctx$nt)) for (c in seq_len(r)) {
      index <- match(paste(factor, r, c, sep = ":"), keys)
      if (is.na(index)) {
        if (!dmu_is_fixed(ctx, role, r, c)) ge_stop("A free DMU component is missing", "incomplete_output")
      } else {
        value <- raw[index, 4L]
        if (dmu_is_fixed(ctx, role, r, c) && abs(value) > 1e-10) {
          ge_stop("DMU did not apply the configured zero covariance", "constraint_mismatch")
        }
        matrix[r, c] <- matrix[c, r] <- value
      }
    }
    if (any(diag(matrix) < 0) || min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) <
        -1e-6 * max(1, max(abs(matrix)))) {
      ge_stop("DMU final covariance matrix is outside the parameter space", "invalid_fit")
    }
    components[[dmu_block_key(role)]] <- matrix
  }
  # R5.6 writes PAROUT in the exact lower-triangular free-parameter order.
  for (i in seq_len(nrow(raw))) {
    role <- roles[raw[i, 1L]]
    metadata <- rbind(metadata, data.frame(role = role, trait1 = row[i], trait2 = column[i],
                                           estimate = raw[i, 4L], factor = raw[i, 1L],
                                           fixed = dmu_is_fixed(ctx, role, row[i], column[i])))
  }
  list(components = components, metadata = metadata)
}

dmu_parse_sampling <- function(ctx, path, metadata) {
  warnings <- character()
  native_se <- setNames(lapply(dmu_roles(ctx), function(role) matrix(NA_real_, ctx$nt, ctx$nt)),
                         vapply(dmu_roles(ctx), dmu_block_key, character(1)))
  free <- !metadata$fixed
  failure <- function(reason) list(metadata = metadata[free, c("role", "trait1", "trait2")],
                                  covariance = NULL, native_se = native_se,
                                  warnings = paste("DMU sampling covariance unavailable:", reason))
  if (!file.exists(path)) return(failure("PAROUT_STD was not produced"))
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(failure("PAROUT_STD is empty"))
  npar <- dmu_number(lines[1L])
  if (length(npar) != 1L || !is.finite(npar) || npar != nrow(metadata) || npar != floor(npar)) {
    return(failure("parameter count does not match the component manifest including restricted zeros"))
  }
  expected <- 1L + npar + npar * (npar + 1L) / 2L
  if (length(lines) != expected) return(failure("PAROUT_STD is truncated or malformed"))
  parse_part <- function(x, width) {
    data <- lapply(strsplit(x, "[[:space:]]+"), dmu_number)
    if (any(lengths(data) != width) || any(!is.finite(unlist(data)))) return(NULL)
    do.call(rbind, data)
  }
  estimates <- parse_part(lines[seq.int(2L, 1L + npar)], 3L)
  pairs <- parse_part(lines[seq.int(2L + npar, expected)], 4L)
  if (is.null(estimates) || is.null(pairs) ||
      !identical(as.integer(estimates[, 1L]), seq_len(npar)) || any(estimates[, 3L] < 0)) {
    return(failure("invalid parameter indices or standard errors"))
  }
  if (any(abs(estimates[, 2L] - metadata$estimate) > 1e-6 * pmax(1, abs(metadata$estimate)))) {
    return(failure("PAROUT_STD order or estimates differ from PAROUT"))
  }
  for (i in seq_len(npar)) {
    key <- dmu_block_key(metadata$role[i])
    r <- metadata$trait1[i]
    c <- metadata$trait2[i]
    native_se[[key]][r, c] <- native_se[[key]][c, r] <- estimates[i, 3L]
  }
  if (any(pairs[, 1:2] != floor(pairs[, 1:2])) ||
      any(pairs[, 1:2] < 1 | pairs[, 1:2] > npar) || any(abs(pairs[, 3L]) > 1 + 1e-6)) {
    return(failure("invalid sampling covariance indices or correlations"))
  }
  keys <- paste(pmax(pairs[, 1L], pairs[, 2L]), pmin(pairs[, 1L], pairs[, 2L]), sep = ":")
  if (anyDuplicated(keys)) return(failure("duplicate parameter pairs"))
  covariance <- matrix(NA_real_, npar, npar)
  for (i in seq_len(nrow(pairs))) {
    r <- pairs[i, 1L]
    c <- pairs[i, 2L]
    covariance[r, c] <- covariance[c, r] <- pairs[i, 4L]
    implied <- pairs[i, 3L] * estimates[r, 3L] * estimates[c, 3L]
    if (abs(implied - pairs[i, 4L]) > 2e-5 * max(1, abs(implied))) {
      return(failure("covariances disagree with their correlations and standard errors"))
    }
  }
  if (anyNA(covariance) || any(diag(covariance) < 0) ||
      any(abs(diag(covariance) - estimates[, 3L]^2) > 1e-5 * pmax(1, estimates[, 3L]^2))) {
    return(failure("invalid diagonal or incomplete sampling matrix"))
  }
  eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -1e-6 * max(1, max(abs(eigenvalues)))) {
    return(failure("sampling covariance is not positive semidefinite"))
  }
  if (any(metadata$fixed) &&
      any(abs(covariance[metadata$fixed, , drop = FALSE]) > 1e-10)) {
    return(failure("restricted parameters have nonzero sampling covariance"))
  }
  result_meta <- metadata[free, c("role", "trait1", "trait2")]
  result_meta$SE <- estimates[free, 3L]
  list(metadata = result_meta, covariance = covariance[free, free, drop = FALSE],
       native_se = native_se, warnings = warnings)
}

dmu_fit <- function(ctx, directory, initial = ctx$initial) {
  stage <- dmu_stage(ctx, directory, initial, task = 1L)
  execution <- dmu_run_backend(ctx, stage, "fit")
  parsed <- dmu_parse_components(ctx, file.path(stage$directory, "PAROUT"))
  sampling <- dmu_parse_sampling(ctx, file.path(stage$directory, "PAROUT_STD"), parsed$metadata)
  converged <- any(grepl("^[[:space:]]*Converged on", execution$output))
  rounds <- grep("^[[:space:]]*[0-9]+[[:space:]]+[-+0-9.EeDd]+.*\\|", execution$output, value = TRUE)
  rounds <- if (length(rounds)) max(vapply(strsplit(trimws(rounds), "[[:space:]]+"),
                                          function(x) as.integer(x[1L]), integer(1))) else NA_integer_
  likelihood <- grep("-2LogL", execution$output, value = TRUE, fixed = TRUE)
  minus2logL <- if (length(likelihood)) {
    value <- sub("^.*?-2LogL[[:space:].]*", "", likelihood[length(likelihood)], perl = TRUE)
    dmu_number(strsplit(trimws(value), "[[:space:]]+")[[1L]][1L])
  } else NA_real_
  warnings <- c(sampling$warnings,
                "DMU reports -2LogL plus a data/fixed-effects constant; raw likelihoods are not compared across engines.")
  if (!converged) warnings <- c(warnings, "DMU did not report numerical convergence.")
  if (any(grepl("Converged on Hessian", execution$output, fixed = TRUE))) {
    warnings <- c(warnings, "DMU accepted its Hessian convergence criterion; update and gradient norms may remain above their thresholds.")
  }
  state <- list(stage = stage$directory, version = execution$version,
                 elapsed = execution$elapsed, task = "AI-REML", inbreeding = isTRUE(ctx$inbreeding),
                 stable_reml = isTRUE(ctx$stable_reml), method = stage$manifest$method,
                 likelihood_constant = TRUE, parameter_covariance_source = "PAROUT_STD")
  ge_write_json(state, file.path(stage$directory, "backend_state.json"))
  list(engine = "DMU", components = parsed$components, param_meta = sampling$metadata,
       param_cov = sampling$covariance, converged = converged, rounds = rounds,
       minus2logL = minus2logL, native_se = sampling$native_se, warnings = warnings,
       files = list.files(stage$directory, full.names = TRUE), version = execution$version,
       state = state)
}

dmu_parse_inbreeding <- function(ctx, directory) {
  if (!isTRUE(ctx$inbreeding)) return(setNames(rep(NA_real_, nrow(ctx$ped)), as.character(ctx$ped$id)))
  path <- file.path(directory, "INBREED")
  if (!file.exists(path) || file.info(path)$size == 0) {
    ge_stop("DMU PED 1 did not produce verified inbreeding coefficients", "missing_inbreeding_output")
  }
  raw <- dmu_numeric_rows(path, 9L)
  id <- as.integer(raw[, 1L])
  if (anyDuplicated(id) || !setequal(id, ctx$ped$id)) {
    ge_stop("DMU inbreeding IDs or coverage do not match the common pedigree", "inbreeding_mismatch")
  }
  # The verified Windows R5.6 build writes the diagonal of A (1+F), despite
  # the manual calling this field F. Founders identify the scale unambiguously.
  founders <- ctx$ped$id[ctx$ped$sire == 0 & ctx$ped$dam == 0]
  founder_values <- raw[match(founders, id), 8L]
  if (!length(founder_values)) ge_stop("DMU F scale cannot be verified without pedigree founders", "unverified_inbreeding")
  if (all(abs(founder_values - 1) < 1e-6)) offset <- 1 else
    if (all(abs(founder_values) < 1e-6)) offset <- 0 else
      ge_stop("DMU INBREED scale is inconsistent with known founders", "unverified_inbreeding")
  values <- raw[, 8L] - offset
  if (any(values < -1e-8 | values > 1 + 1e-8)) ge_stop("Invalid DMU inbreeding coefficients", "inbreeding_mismatch")
  setNames(values, as.character(id))
}

dmu_parse_predictions <- function(ctx, directory) {
  path <- file.path(directory, "RESIDUAL")
  if (!file.exists(path)) return(NULL)
  # Verified R5.6 direct-solver output: record, observed, predicted, residual,
  # hat-matrix diagonal, observed corrected for fixed effects (nt of each).
  raw <- dmu_numeric_rows(path, 1L + 5L * ctx$nt)
  row <- as.integer(raw[, 1L])
  if (nrow(raw) != nrow(ctx$data) || anyDuplicated(row) ||
      !setequal(row, seq_len(nrow(ctx$data)))) {
    ge_stop("DMU prediction record indices do not match the staged observations", "invalid_output")
  }
  output <- lapply(seq_len(ctx$nt), function(trait) {
    predicted <- raw[, 1L + ctx$nt + trait]
    observed <- ctx$data[[paste0("y", trait)]][row]
    present <- !is.na(observed)
    native_observed <- raw[, 1L + trait]
    native_residual <- raw[, 1L + 2L * ctx$nt + trait]
    tolerance <- 2e-6 * pmax(1, abs(observed[present]))
    if (any(abs(native_observed[present] - observed[present]) > tolerance) ||
        any(abs(predicted[present] + native_residual[present] - observed[present]) > tolerance)) {
      ge_stop("DMU observed/predicted/residual columns fail the original-observation check", "invalid_output")
    }
    predicted[!present] <- NA_real_
    data.frame(record_id = ctx$data$record_id[row], trait = trait, predicted = predicted)
  })
  do.call(rbind, output)
}

dmu_blup <- function(ctx, fit, directory, components = fit$components) {
  stage <- dmu_stage(ctx, directory, components, task = 11L)
  execution <- dmu_run_backend(ctx, stage, "blup")
  path <- file.path(stage$directory, "SOL")
  if (!file.exists(path)) ge_stop("DMU4 did not produce SOL", "missing_output")
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  header <- dmu_number(strsplit(lines[1L], "[[:space:]]+")[[1L]])
  if (length(header) != 9L || header[1L] != 0L || header[3L] != 7L || header[4L] != 2L) {
    ge_stop("Unexpected DMU4 direct-solver SOL descriptor; SEP cannot be inferred", "unsupported_output_format")
  }
  raw <- dmu_numeric_rows(path, 9L)[-1L, , drop = FALSE]
  animal <- raw[raw[, 1L] == 4L & raw[, 4L] == 1L, , drop = FALSE]
  if (!nrow(animal) || any(animal[, 2L] != animal[, 3L]) ||
      any(animal[, 2L] < 1L | animal[, 2L] > ctx$nt)) {
    ge_stop("DMU additive-factor trait mapping is invalid", "invalid_output")
  }
  ids <- as.integer(animal[, 5L])
  traits <- as.integer(animal[, 2L])
  expected <- expand.grid(id = ctx$ped$id, trait = seq_len(ctx$nt))
  keys <- paste(ids, traits, sep = ":")
  expected_keys <- paste(expected$id, expected$trait, sep = ":")
  if (anyDuplicated(keys) || !setequal(keys, expected_keys)) {
    ge_stop(paste("DMU additive EBV coverage is incomplete or unexpected:", length(keys),
                  "rows found,", length(expected_keys), "expected. Check version and authorization."),
            "incomplete_ebv_output")
  }
  if (any(animal[, 9L] < 0)) ge_stop("DMU SEP must be nonnegative", "invalid_output")
  inbreeding <- dmu_parse_inbreeding(ctx, stage$directory)
  # Build one observed-animal set per trait; ancestors and missing traits remain FALSE.
  has_record <- logical(length(ids))
  for (trait in seq_len(ctx$nt)) {
    present_ids <- unique(ctx$data$id[!is.na(ctx$data[[paste0("y", trait)]])])
    rows <- traits == trait
    has_record[rows] <- ids[rows] %in% present_ids
  }
  result <- data.frame(id = ids, original_id = ctx$ped$original_id[match(ids, ctx$ped$id)],
                        trait = traits, EBV = animal[, 8L], SEP = animal[, 9L],
                        F = unname(inbreeding[as.character(ids)]), has_record = has_record,
                        stringsAsFactors = FALSE)
  if (anyNA(result$original_id) || (isTRUE(ctx$inbreeding) && anyNA(result$F))) {
    ge_stop("DMU original-ID or verified-inbreeding mapping is incomplete", "invalid_output")
  }
  result <- result[order(result$id, result$trait), , drop = FALSE]
  rownames(result) <- NULL
  ge_write_json(list(version = execution$version, task = "DMU4 FSPAK method 9",
                      predicted_rows = nrow(result), SEP_source = "SOL second real field",
                      inbreeding_source = if (isTRUE(ctx$inbreeding)) "INBREED" else "PED 2: inbreeding not estimated or verified; F=NA"),
                file.path(stage$directory, "backend_state.json"))
  attr(result, "engine_files") <- list.files(stage$directory, full.names = TRUE)
  attr(result, "fitted") <- dmu_parse_predictions(ctx, stage$directory)
  attr(result, "warnings") <- if (isTRUE(ctx$inbreeding))
    "DMU R5.6 INBREED values are rounded in the public text file; F is verified from founder scaling and may have approximately 0.0005 rounding error." else
    "Inbreeding is disabled (PED 2); F is unverified and returned as NA. Reliability requires verified F and SEP."
  result
}
