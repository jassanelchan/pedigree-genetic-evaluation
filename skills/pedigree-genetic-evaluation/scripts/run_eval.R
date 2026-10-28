arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
script_directory <- dirname(script_path)
skill_directory <- dirname(script_directory)

for (module in c("common.R", "data.R", "f90.R", "dmu.R", "results.R")) {
  file <- file.path(script_directory, module)
  if (file.exists(file)) source(file, encoding = "UTF-8")
}

ge_option <- function(name, default = NULL) {
  index <- match(name, arguments)
  if (is.na(index)) return(default)
  if (index >= length(arguments)) ge_stop(paste("Missing value after", name))
  arguments[[index + 1L]]
}

ge_validate_spec <- function(spec) {
  check <- function(value, allowed, label) {
    unknown <- setdiff(names(value), allowed)
    if (length(unknown)) ge_stop(paste(label, "contains unsupported fields:", paste(unknown,collapse=", ")),
                                "unsupported_spec_field")
  }
  source_check <- function(value, allow_supplements=FALSE) {
    if (is.null(value)) return(invisible(NULL))
    check(value,c("path","sheet","columns","missing_values","date_origin","date_format","header","sep","encoding","supplements"),"Table source")
    check(value$columns,c("id","sire","dam"),"ID/parent column map")
    if (length(value$supplements) && !allow_supplements) ge_stop("Parent supplements belong only to the main pedigree source; nested/phenotype supplements are unsupported")
    for (extra in value$supplements) source_check(extra)
  }
  check(spec,c("schema_version","engine","output_dir","software","timeout","threads","pedigree","analyses"),"Project")
  check(spec$software,c("renumf90","blupf90","dmu_dir"),"Software")
  source_check(spec$pedigree,TRUE)
  for (analysis in spec$analyses) {
    check(analysis,c("name","species","model","phenotypes","pedigree","traits","fixed","shared","event_date","birth_date",
                     "record_key","qc","inbreeding","predict_ids","allow_observed_founders","trend_basis","initial","zero_covariances","random"),"Analysis")
    source_check(analysis$phenotypes)
    source_check(analysis$pedigree,TRUE)
    for (trait in analysis$traits) check(trait,c("name","label","column","unit","distribution","record_type","record_filter"),"Trait")
    for (effect in analysis$fixed) check(effect,c("name","column","type","from_date","grouping","traits"),"Fixed effect")
    for (trait in analysis$traits) {
      check(trait$record_filter, c("column", "values"), "Trait record filter")
    }
    for (term in analysis$random) {
      check(term, c("name", "role", "column", "dam", "traits"), "Random effect")
    }
    check(analysis$shared,c("column","dam"),"Shared effect")
    check(analysis$qc,c("method","sd","group_by"),"QC")
    for (zero in analysis$zero_covariances) check(zero,c("component","trait1","trait2"),"Zero covariance")
    if (!is.null(analysis$trend_basis) && !analysis$trend_basis %in% c("event_year","birth_year")) {
      ge_stop("trend_basis must be event_year or birth_year")
    }
    for (effect in analysis$fixed) if (!is.null(effect$from_date) && !effect$from_date %in% c("event","birth")) {
      ge_stop("from_date must be event or birth")
    }
  }
  if (!is.null(spec$timeout) && (!is.finite(spec$timeout) || spec$timeout <= 0)) ge_stop("Timeout must be positive seconds")
  if (!is.null(spec$threads) && (!is.finite(spec$threads) || spec$threads < 1L || spec$threads != floor(spec$threads))) ge_stop("Threads must be a positive integer")
  invisible(TRUE)
}

ge_snapshot_inputs <- function(spec, base, directory, software) {
  snapshots <- file.path(directory,"input_sources")
  dir.create(snapshots,recursive=TRUE,showWarnings=FALSE)
  rows <- list()
  seen <- new.env(parent=emptyenv())
  snapshot <- function(source) {
    if (is.null(source)) return(NULL)
    original <- ge_path(source$path,base)
    hash <- ge_hash_file(original)
    name <- if(exists(hash,envir=seen,inherits=FALSE)) get(hash,envir=seen) else paste0(substr(hash,1L,20L),"_",basename(original))
    assign(hash,name,envir=seen)
    destination <- file.path(snapshots,name)
    if (!file.exists(destination) && !file.copy(original,destination,overwrite=FALSE)) {
      ge_stop(paste("Cannot preserve source input:",original),"input_snapshot_failed")
    }
    if (!identical(ge_hash_file(destination),hash)) ge_stop("Input snapshot SHA256 mismatch","input_snapshot_failed")
    rows[[length(rows)+1L]] <<- data.frame(original_path=original,snapshot_path=paste0("input_sources/",name),
                                         sha256=hash,bytes=file.info(original)$size)
    source$path <- paste0("input_sources/",name)
    if (length(source$supplements)) source$supplements <- lapply(source$supplements,snapshot)
    source
  }
  reproduced <- spec
  reproduced$software <- software
  reproduced$pedigree <- snapshot(spec$pedigree)
  for (k in seq_along(spec$analyses)) {
    reproduced$analyses[[k]]$phenotypes <- snapshot(spec$analyses[[k]]$phenotypes)
    reproduced$analyses[[k]]$pedigree <- snapshot(spec$analyses[[k]]$pedigree)
  }
  reproduced$output_dir <- "reproduced_runs"
  ge_write_json(reproduced,file.path(directory,"reproduce_project.json"))
  manifest <- unique(do.call(rbind,rows))
  ge_write_bom_csv(manifest,file.path(directory,"input_manifest.csv"))
  manifest
}

ge_software_profile <- function(spec) {
  profile_path <- file.path(skill_directory, "host_profile.json")
  result <- if (file.exists(profile_path)) ge_read_json(profile_path)$software else list()
  for (name in names(spec$software)) result[[name]] <- spec$software[[name]]
  for (name in c("blupf90", "renumf90")) {
    if (is.null(result[[name]])) {
      detected <- Sys.which(if (name == "blupf90") "blupf90+" else "renumf90")
      if (nzchar(detected)) result[[name]] <- unname(detected)
    }
  }
  if (!is.null(result$dmu_dir) && !file.exists(file.path(result$dmu_dir,"dmu1.exe")) &&
      file.exists(file.path(result$dmu_dir,"bin","dmu1.exe"))) result$dmu_dir <- file.path(result$dmu_dir,"bin")
  result
}

ge_fit_validity <- function(ctx, fit) {
  if (!isTRUE(fit$converged)) ge_stop("REML did not reach a supported convergence condition", "not_converged")
  roles <- c(vapply(ge_random_terms(ctx), function(term) term$name, character(1)), "residual")
  for (role in roles) {
    matrix <- fit$components[[role]]
    if (is.null(matrix) || !all(dim(matrix) == c(ctx$nt, ctx$nt)) || any(!is.finite(matrix))) {
      ge_stop(paste("Incomplete final component:", role), "incomplete_components")
    }
    scale <- max(1, max(abs(matrix)))
    if (!is.null(ctx$constraints[[role]]) &&
        any(abs(matrix[ctx$constraints[[role]]]) > 1e-8 * scale)) {
      ge_stop(paste("Fixed or structural zero changed in final component:", role), "invalid_covariance")
    }
    if (max(abs(matrix - t(matrix))) > 1e-8 * scale ||
        min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * scale) {
      ge_stop(paste("Invalid final covariance matrix:", role), "invalid_covariance")
    }
  }
  invisible(TRUE)
}

ge_execute_fit <- function(ctx, engine, directory, resume) {
  backend_file <- file.path(script_directory, if (engine == "blupf90") "f90.R" else "dmu.R")
  executable_paths <- if (engine == "blupf90") c(ctx$software$blupf90, ctx$software$renumf90) else {
    if (is.null(ctx$software$dmu_dir)) character() else file.path(ctx$software$dmu_dir,
                                                               c("dmu1.exe", "dmuai.exe", "dmu4.exe"))
  }
  expected_programs <- if (engine == "blupf90") 2L else 3L
  if (length(executable_paths) != expected_programs || any(!file.exists(executable_paths))) {
    ge_stop(paste("Executable environment is unavailable for", engine), "backend_unavailable")
  }
  observation_columns <- c("record_id", "id", paste0("y", seq_len(ctx$nt)),
                            vapply(ctx$fixed, function(x) x$column, character(1)),
                            vapply(ge_random_terms(ctx), function(term) term$column, character(1)))
  observation_columns <- unique(observation_columns)
  fingerprint <- ge_hash(list(data = ctx$data[observation_columns], ped = ctx$ped,
                              fixed = ctx$fixed, random = ge_random_terms(ctx), constraints = ctx$constraints,
                              inbreeding = ctx$inbreeding,
                              initial = ctx$initial, family = ctx$family,
                              programs = vapply(executable_paths, ge_hash_file, character(1)),
                              backend = ge_hash_file(backend_file),
                              common = ge_hash_file(file.path(script_directory, "common.R")),
                              preprocess = ge_hash_file(file.path(script_directory, "data.R"))))
  cache_path <- file.path(directory, "verified_fit.rds")
  if (resume && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$fingerprint, fingerprint)) {
      ge_fit_validity(ctx, cached$fit)
      cat("Reusing verified fit:", ctx$name, engine, "\n")
      return(cached)
    }
  }
  attempts <- list()
  fit_function <- if (engine == "blupf90") f90_fit else dmu_fit
  fit <- NULL
  for (attempt in seq_len(3L)) {
    initial <- ctx$initial
    if (attempt == 2L) {
      initial$additive <- initial$additive * (0.5 / 0.3)
      initial$residual <- initial$residual * 0.8
      signs <- (-1) ^ seq_len(ctx$nt)
      for (role in names(initial)) {
        initial[[role]] <- initial[[role]] * outer(signs, signs)
      }
    }
    context <- ctx
    context$attempt <- attempt
    context$stable_reml <- attempt == 3L
    attempt_path <- ge_stage_directory(directory, paste0("reml_attempt_", attempt))
    cat("Running REML:", ctx$name, engine, "attempt", attempt, "\n")
    outcome <- tryCatch({
      value <- fit_function(context, attempt_path, initial)
      ge_fit_validity(context, value)
      value
    }, error = function(e) e)
    if (!inherits(outcome, "error")) {
      fit <- outcome
      attempts[[length(attempts) + 1L]] <- list(attempt = attempt, directory = attempt_path,
                                               state = "converged", warnings = fit$warnings)
      break
    }
    attempts[[length(attempts) + 1L]] <- list(attempt = attempt, directory = attempt_path,
                                             state = "failed", reason = conditionMessage(outcome),
                                             code = if (inherits(outcome, "ge_error")) outcome$code else "runtime_error")
    if (inherits(outcome, "ge_error") && outcome$code %in%
        c("backend_unavailable", "unsupported_version", "unsupported_engine_version", "invalid_input", "unsupported_constraint", "incomplete_parameter_map")) break
  }
  ge_write_json(attempts, file.path(directory, "attempts.json"))
  if (is.null(fit)) ge_stop(paste("REML failed after bounded attempts:", attempts[[length(attempts)]]$reason),
                          "fit_failed")
  result <- list(fingerprint = fingerprint, fit = fit, attempts = attempts)
  saveRDS(result, cache_path)
  result
}

ge_ebv_coverage <- function(ctx, ebv) {
  columns <- c("id", "original_id", "trait", "EBV", "SEP", "F", "has_record")
  if (!all(columns %in% names(ebv))) ge_stop("BLUP result lacks required columns", "incomplete_ebv")
  if (nrow(ebv) != nrow(ctx$ped) * ctx$nt ||
      anyDuplicated(ebv[c("id", "trait")]) || any(!is.finite(ebv$EBV)) ||
      !setequal(ebv$id, ctx$ped$id) || !setequal(ebv$trait, seq_len(ctx$nt))) {
    ge_stop("Target animal by trait output coverage is incomplete", "incomplete_ebv")
  }
  index <- match(ebv$id, ctx$ped$id)
  if (anyNA(index) || any(ebv$original_id != ctx$ped$original_id[index])) {
    ge_stop("BLUP result ID mapping mismatch", "invalid_id_map")
  }
  invisible(TRUE)
}

ge_compare_results <- function(ctx, f90_result, dmu_result, directory) {
  printed_bound <- function(value, digits) {
    result <- rep(0, length(value))
    valid <- is.finite(value) & value != 0
    result[valid] <- 0.5 * 10^(floor(log10(abs(value[valid]))) - digits + 1L)
    result
  }
  rows <- list()
  for (role in names(f90_result$fit$components)) {
    first <- f90_result$fit$components[[role]]
    second <- dmu_result$fit$components[[role]]
    for (i in seq_len(ctx$nt)) for (j in seq_len(ctx$nt)) {
      tolerance <- 1e-6 + 1e-4 * max(abs(first[i,j]), abs(second[i,j]))
      rows[[length(rows) + 1L]] <- data.frame(component = role, trait1 = i, trait2 = j,
                                             blupf90 = first[i,j], dmu = second[i,j],
                                             absolute_difference = abs(first[i,j] - second[i,j]),
                                             tolerance = tolerance, pass = abs(first[i,j] - second[i,j]) <= tolerance)
    }
  }
  components <- do.call(rbind, rows)
  fixed_path <- ge_stage_directory(directory, "fixed_component_dmu_blup")
  shared <- dmu_blup(ctx, dmu_result$fit, fixed_path, components = f90_result$fit$components)
  ge_ebv_coverage(ctx, shared)
  first <- f90_result$ebv[order(f90_result$ebv$id, f90_result$ebv$trait), ]
  second <- shared[order(shared$id, shared$trait), ]
  rows <- lapply(seq_len(ctx$nt), function(j) {
    mask <- first$trait == j
    difference <- first$EBV[mask] - second$EBV[mask]
    scale <- max(1, sd(ctx$data[[paste0("y",j)]], na.rm = TRUE))
    tolerance <- printed_bound(second$EBV[mask], 6L) + 5e-9 + 2e-7 * scale
    pev_difference <- first$SEP[mask]^2 - second$SEP[mask]^2
    valid_pev <- is.finite(pev_difference)
    sep <- second$SEP[mask]
    sep_bound <- printed_bound(sep, 6L) + 5e-9
    pev_tolerance <- (2 * abs(sep) + sep_bound) * sep_bound +
      2e-7 * max(1, f90_result$fit$components$additive[j,j])
    data.frame(trait = ctx$traits$name[j], n = sum(mask), max_absolute_EBV_difference = max(abs(difference)),
                mean_EBV_difference = mean(difference),
                EBV_correlation = if (sd(first$EBV[mask]) > 0 && sd(second$EBV[mask]) > 0) cor(first$EBV[mask],second$EBV[mask]) else NA_real_,
                EBV_tolerance = max(tolerance), EBV_pass = all(abs(difference) <= tolerance),
                max_absolute_PEV_difference = if (any(valid_pev)) max(abs(pev_difference[valid_pev])) else NA_real_,
                PEV_coverage = sum(valid_pev), PEV_tolerance = max(pev_tolerance),
                PEV_pass = all(valid_pev) && all(abs(pev_difference) <= pev_tolerance),
                stringsAsFactors = FALSE)
  })
  predictions <- do.call(rbind, rows)
  fitted_first <- attr(f90_result$ebv, "fitted")
  fitted_second <- attr(shared, "fitted")
  fitted_comparison <- NULL
  if (!is.null(fitted_first) && !is.null(fitted_second)) {
    joined <- merge(fitted_first, fitted_second, by = c("record_id", "trait"), suffixes = c("_f90", "_dmu"))
    joined$difference <- joined$predicted_f90 - joined$predicted_dmu
    joined$tolerance <- printed_bound(joined$predicted_f90, 8L) +
      2e-7 * pmax(1, abs(joined$predicted_f90))
    joined$pass <- is.finite(joined$difference) & abs(joined$difference) <= joined$tolerance
    observed <- sum(vapply(seq_len(ctx$nt), function(j) sum(!is.na(ctx$data[[paste0("y",j)]])), integer(1)))
    valid <- is.finite(joined$predicted_f90) & is.finite(joined$predicted_dmu)
    fitted_pass <- sum(valid) == observed && all(joined$pass[valid])
    fitted_comparison <- joined
  } else fitted_pass <- NA
  ge_write_bom_csv(components, file.path(directory, "component_comparison.csv"))
  ge_write_bom_csv(predictions, file.path(directory, "fixed_component_EBV_comparison.csv"))
  if (!is.null(fitted_comparison)) ge_write_bom_csv(fitted_comparison, file.path(directory, "fitted_comparison.csv"))
  list(components = components, predictions = predictions,
        fitted = fitted_comparison, fitted_verified = fitted_pass,
        precision_policy = "DMU SOL 6 significant digits; F90 solutions 8 decimals; F90 fitted 8 significant digits. Propagate print bounds plus a small numerical allowance for native input precision.",
        pass = all(components$pass) && all(predictions$EBV_pass) && all(predictions$PEV_pass) && !identical(fitted_pass, FALSE))
}

ge_main <- function() {
  ge_require(c("jsonlite", "digest", "processx", "openxlsx", "Matrix"))
  spec_file <- ge_option("--spec")
  if (is.null(spec_file)) ge_stop("Usage: Rscript --vanilla run_eval.R --spec project.json [--output DIR] [--resume RUN_DIR] [--prepare-only]",
                                "missing_project_spec")
  spec_file <- ge_path(spec_file)
  base <- dirname(spec_file)
  spec <- ge_read_json(spec_file)
  ge_validate_spec(spec)
  schema <- suppressWarnings(as.numeric(unlist(spec$schema_version)))
  if (length(schema) != 1L || !is.finite(schema) || !schema %in% c(1L, 2L)) {
    ge_stop("schema_version must be integer 1 or 2", "unsupported_schema")
  }
  engine <- ge_option("--engine", spec$engine)
  if (is.null(engine) || !engine %in% c("blupf90", "dmu", "compare")) {
    ge_stop("Choose blupf90, dmu or compare explicitly", "engine_choice_required")
  }
  if (!length(spec$analyses)) ge_stop("At least one analysis is required")
  names <- vapply(spec$analyses, function(x) if (is.null(x$name)) "" else x$name, character(1))
  if (any(!nzchar(names)) || anyDuplicated(names)) ge_stop("Analysis names must be nonempty and unique")
  if (any(vapply(spec$analyses, function(x) is.null(x$species) || !nzchar(x$species), logical(1)))) {
    ge_stop("Each analysis needs a known species; do not pool species", "species_required")
  }
  output <- ge_option("--output", spec$output_dir)
  if (is.null(output)) output <- file.path(base, "genetic_evaluation_runs")
  output <- ge_path(output, base, must_work = FALSE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  resume_directory <- ge_option("--resume")
  resume <- !is.null(resume_directory)
  run_directory <- if (resume) ge_path(resume_directory) else {
    ge_stage_directory(output, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", substr(ge_hash(spec),1L,8L)))
  }
  assign("ge_current_run_directory", run_directory, envir = .GlobalEnv)
  delivery <- file.path(run_directory, "delivery")
  dir.create(delivery, recursive = TRUE, showWarnings = FALSE)
  code_directory <- ge_stage_directory(file.path(run_directory, "code_records"),
                                       paste0("invocation_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  for (file in list.files(script_directory, pattern = "\\.R$", full.names = TRUE)) {
    file.copy(file, file.path(code_directory, basename(file)), overwrite = TRUE)
  }
  file.copy(spec_file, file.path(code_directory, "project_spec.json"), overwrite = TRUE)
  software <- ge_software_profile(spec)
  input_snapshots <- ge_snapshot_inputs(spec,base,run_directory,software)
  manifest <- list(schema_version = 1L, engine = engine, state = "running",
                    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                    spec_path = spec_file, spec_sha256 = ge_hash_file(spec_file),
                    software = software, R_version = as.character(getRversion()), analyses = list())
  manifest$input_snapshots <- input_snapshots
  manifest_file <- file.path(run_directory, "run_manifest.json")
  ge_write_json(manifest, manifest_file)
  workbook <- openxlsx::createWorkbook()
  errors <- list()
  compare_failed <- FALSE
  parameters <- list()
  descriptions <- list()
  qc_tables <- list()
  coverage_tables <- list()
  input_tables <- list()
  model_rows <- list()
  annual_rows <- list()
  trend_rows <- list()
  contexts <- list()
  engines <- if (engine == "compare") c("blupf90", "dmu") else engine
  prepare_only <- "--prepare-only" %in% arguments
  for (index in seq_along(spec$analyses)) {
    analysis <- spec$analyses[[index]]
    analysis_directory <- file.path(run_directory, paste0("analysis_", index))
    dir.create(analysis_directory, recursive = TRUE, showWarnings = FALSE)
    cat("Preparing:", analysis$name, "\n")
    context <- tryCatch(ge_prepare_context(spec, analysis, base, software), error = function(e) e)
    if (inherits(context, "error")) {
      errors[[length(errors) + 1L]] <- list(analysis = analysis$name, phase = "prepare", reason = conditionMessage(context),
                                          code = if (inherits(context,"ge_error")) context$code else "runtime_error")
      next
    }
    contexts[[analysis$name]] <- context
    saveRDS(context, file.path(analysis_directory, "prepared_context.rds"))
    ge_write_bom_csv(context$ped, file.path(analysis_directory, "animal_id_map.csv"))
    ge_write_bom_csv(context$qc_audit, file.path(analysis_directory, "qc_audit.csv"))
    ge_write_bom_csv(context$fixed_coverage, file.path(analysis_directory, "fixed_effect_coverage.csv"))
    ge_write_bom_csv(context$input_audit, file.path(analysis_directory, "input_audit.csv"))
    ge_write_bom_csv(context$rank_checks, file.path(analysis_directory, "design_rank.csv"))
    pedigree_source <- if (is.null(analysis$pedigree)) spec$pedigree else analysis$pedigree
    supplemental_sources <- lapply(pedigree_source$supplements, function(x) {
      list(path=ge_path(x$path,base),sha256=ge_hash_file(ge_path(x$path,base)))
    })
    ge_write_json(list(model = context$model_description, inbreeding = context$inbreeding,
                        traits = context$traits, fixed = context$fixed, random = context$random,
                        constraints = context$constraints,
                        pedigree_supplements = supplemental_sources,
                        phenotype_sha256 = ge_hash_file(ge_path(analysis$phenotypes$path, base)),
                        pedigree_sha256 = ge_hash_file(ge_path(if (is.null(analysis$pedigree)) spec$pedigree$path else analysis$pedigree$path, base))),
                   file.path(analysis_directory, "analysis_definition.json"))
    descriptions[[length(descriptions)+1L]] <- cbind(analysis = analysis$name, context$descriptive)
    qc_tables[[length(qc_tables)+1L]] <- cbind(analysis = analysis$name, context$qc_audit)
    coverage_tables[[length(coverage_tables)+1L]] <- cbind(analysis = analysis$name, context$fixed_coverage)
    input_tables[[length(input_tables)+1L]] <- cbind(analysis = analysis$name, context$input_audit)
    model_rows[[length(model_rows)+1L]] <- data.frame(analysis = analysis$name, species = context$species,
                                                     model = context$model_description, inbreeding = context$inbreeding,
                                                     n_records = nrow(context$data), n_animals = length(unique(context$data$id)),
                                                     pedigree_animals = nrow(context$ped), trend_basis = context$trend_basis)
    if (prepare_only) next
    finished <- list()
    for (backend in engines) {
      backend_directory <- file.path(analysis_directory, backend)
      dir.create(backend_directory, recursive = TRUE, showWarnings = FALSE)
      fit <- NULL
      phase <- "fit"
      outcome <- tryCatch({
        cached <- ge_execute_fit(context, backend, backend_directory, resume)
        fit <- cached$fit
        phase <- "blup"
        cache_ebv <- file.path(backend_directory, "verified_ebv.rds")
        if (resume && file.exists(cache_ebv)) {
          saved <- readRDS(cache_ebv)
          ebv <- if (identical(saved$fingerprint, cached$fingerprint)) saved$ebv else NULL
        } else ebv <- NULL
        if (is.null(ebv)) {
          prediction_directory <- ge_stage_directory(backend_directory, "blup_final")
          cat("Predicting EBV:", analysis$name, backend, "\n")
          function_blup <- if (backend == "blupf90") f90_blup else dmu_blup
          ebv <- function_blup(context, fit, prediction_directory)
          ge_ebv_coverage(context, ebv)
          saveRDS(list(fingerprint = cached$fingerprint, ebv = ebv), cache_ebv)
        } else ge_ebv_coverage(context, ebv)
        result_directory <- file.path(delivery, paste0("analysis_", index, "_", backend))
        dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
        phase <- "export"
        exported <- ge_export_results(context, fit, ebv, result_directory)
        values <- exported$parameters
        values$analysis <- analysis$name
        values$engine <- backend
        parameters[[length(parameters)+1L]] <- values
        if (!is.null(exported$annual) && nrow(exported$annual)) annual_rows[[length(annual_rows)+1L]] <- cbind(analysis=analysis$name,engine=backend,exported$annual)
        if (!is.null(exported$trends) && nrow(exported$trends)) trend_rows[[length(trend_rows)+1L]] <- cbind(analysis=analysis$name,engine=backend,exported$trends)
        manifest$analyses[[length(manifest$analyses)+1L]] <- list(name=analysis$name,engine=backend,
                                                                state=if (is.null(fit$param_cov) || anyNA(ebv$SEP)) "partial" else "complete",
                                                                warnings=unique(c(fit$warnings,exported$warnings)),
                                                                fit_fingerprint=cached$fingerprint,
                                                                artifacts=list(exported$workbook,exported$ebv_csv,exported$histograms,exported$trend_figures))
        list(fit=fit,ebv=ebv,fingerprint=cached$fingerprint)
      }, error = function(e) e)
      if (inherits(outcome,"error")) {
        errors[[length(errors)+1L]] <- list(analysis=analysis$name,engine=backend,phase=phase,reason=conditionMessage(outcome),
                                            code=if(inherits(outcome,"ge_error")) outcome$code else "runtime_error")
        if (!is.null(fit) && phase == "blup") {
          values <- ge_parameters(context,fit)
          values$analysis <- analysis$name
          values$engine <- backend
          parameters[[length(parameters)+1L]] <- values
          manifest$analyses[[length(manifest$analyses)+1L]] <- list(name=analysis$name,engine=backend,state="partial",
                                                                  warnings=paste("Verified REML retained; BLUP failed:",conditionMessage(outcome)))
        }
      } else finished[[backend]] <- outcome
      ge_write_json(manifest,manifest_file)
    }
    if (engine == "compare" && all(c("blupf90","dmu") %in% names(finished))) {
      comparison_directory <- file.path(delivery,paste0("analysis_",index,"_comparison"))
      dir.create(comparison_directory,recursive=TRUE,showWarnings=FALSE)
      comparison_fingerprint <- ge_hash(list(f90=finished$blupf90$fingerprint,
                                              dmu=finished$dmu$fingerprint,
                                              traits=context$traits,
                                              comparison_algorithm=body(ge_compare_results)))
      comparison_cache <- file.path(comparison_directory,"verified_comparison.rds")
      compared <- tryCatch({
        value <- NULL
        if (resume && file.exists(comparison_cache)) {
          saved <- readRDS(comparison_cache)
          if (identical(saved$fingerprint,comparison_fingerprint)) value <- saved$result
        }
        if (is.null(value)) {
          value <- ge_compare_results(context,finished$blupf90,finished$dmu,comparison_directory)
          saveRDS(list(fingerprint=comparison_fingerprint,result=value),comparison_cache)
        } else cat("Reusing verified comparison:",analysis$name,"\n")
        value
      },error=function(e)e)
      if (inherits(compared,"error")) {
        errors[[length(errors)+1L]] <- list(analysis=analysis$name,phase="comparison",reason=conditionMessage(compared))
      } else {
        ge_write_json(compared,file.path(comparison_directory,"comparison.json"))
        if (!compared$pass) compare_failed <- TRUE
      }
    }
  }
  add_sheet <- function(name, rows, decimals = 3L) {
    if (!length(rows)) return(invisible(NULL))
    value <- if (is.data.frame(rows)) rows else do.call(rbind,rows)
    openxlsx::addWorksheet(workbook,name)
    openxlsx::writeDataTable(workbook,name,value)
    numeric_columns <- which(vapply(value,is.numeric,logical(1)))
    if (nrow(value) && length(numeric_columns)) {
      style <- openxlsx::createStyle(numFmt=if(decimals == 0L) "0" else paste0("0.",paste(rep("0",decimals),collapse="")))
      openxlsx::addStyle(workbook,name,style,rows=2:(nrow(value)+1L),cols=numeric_columns,gridExpand=TRUE)
      counts <- which(grepl("^n$|^n_|^removed$|^pedigree_animals$|^year$",names(value)))
      if(length(counts)) openxlsx::addStyle(workbook,name,openxlsx::createStyle(numFmt="0"),rows=2:(nrow(value)+1L),cols=counts,gridExpand=TRUE,stack=FALSE)
    }
    openxlsx::setColWidths(workbook,name,cols=seq_len(ncol(value)),widths="auto")
    openxlsx::freezePane(workbook,name,firstRow=TRUE)
  }
  add_sheet("Models",model_rows)
  add_sheet("Descriptive",descriptions,2L)
  add_sheet("QC",qc_tables,2L)
  add_sheet("Fixed_Coverage",coverage_tables,0L)
  add_sheet("Input_Audit",input_tables,0L)
  add_sheet("Variance_Parameters",parameters,3L)
  add_sheet("Annual",annual_rows,3L)
  add_sheet("Trends",trend_rows,3L)
  if(length(errors)) add_sheet("Errors",lapply(errors,function(x)data.frame(analysis=x$analysis,phase=x$phase,
                                                                           engine=if(is.null(x$engine)) "" else x$engine,
                                                                           code=if(is.null(x$code)) "" else x$code,reason=x$reason)))
  if (!length(names(workbook))) {
    openxlsx::addWorksheet(workbook,"Status")
    openxlsx::writeData(workbook,"Status",data.frame(state="No completed analysis; inspect run_manifest and errors"))
  }
  openxlsx::saveWorkbook(workbook,file.path(delivery,"evaluation_summary.xlsx"),overwrite=TRUE)
  manifest$errors <- errors
  manifest$comparison_failed <- compare_failed
  partial <- any(vapply(manifest$analyses,function(x) identical(x$state,"partial"),logical(1)))
  manifest$state <- if(length(errors)||compare_failed||partial) "partial" else if(prepare_only) "prepared" else "complete"
  manifest$finished_utc <- format(Sys.time(),tz="UTC",usetz=TRUE)
  manifest$delivery_directory <- delivery
  manifest$packages <- as.list(vapply(c("jsonlite","digest","processx","openxlsx","Matrix"),function(x)as.character(packageVersion(x)),character(1)))
  ge_write_json(manifest,manifest_file)
  cat("Run state:",manifest$state,"\nOutput:",run_directory,"\n")
  if(length(errors)) for(error in errors) cat("ERROR:",error$analysis,error$reason,"\n")
  if(manifest$state == "partial") 2L else 0L
}

status <- tryCatch(ge_main(), error = function(e) {
  code <- if(inherits(e,"ge_error")) e$code else "runtime_error"
  if(exists("ge_current_run_directory",envir=.GlobalEnv)) {
    ge_write_json(list(state="failed",code=code,message=conditionMessage(e)),
                  file.path(get("ge_current_run_directory",envir=.GlobalEnv),"fatal_error.json"))
  }
  cat("ERROR [",code,"]: ",conditionMessage(e),"\n",sep="",file=stderr())
  2L
})
quit(status=status,save="no")
