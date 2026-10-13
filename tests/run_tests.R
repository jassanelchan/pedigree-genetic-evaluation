arguments <- commandArgs(trailingOnly = TRUE)

script_argument <- commandArgs()[grepl("^--file=", commandArgs())]
script_file <- sub("^--file=", "", script_argument[1L])
repository <- dirname(dirname(normalizePath(script_file, winslash = "/")))
skill_path <- file.path(repository, "skills", "pedigree-genetic-evaluation")

if (length(arguments)) {
  if (length(arguments) != 2L || arguments[1L] != "--skill-path") {
    stop("Usage: Rscript tests/run_tests.R [--skill-path PATH]")
  }

  skill_path <- normalizePath(arguments[2L], winslash = "/", mustWork = TRUE)
}

for (module in c("common.R", "data.R", "results.R")) {
  source(file.path(skill_path, "scripts", module), encoding = "UTF-8")
}

checks <- 0L

check <- function(condition, description) {
  if (!isTRUE(condition)) {
    stop("Failed: ", description, call. = FALSE)
  }

  checks <<- checks + 1L
}

check_close <- function(actual, expected, description, tolerance = 1e-8) {
  check(isTRUE(all.equal(actual, expected, tolerance = tolerance,
                         check.attributes = FALSE)), description)
}

expect_error_code <- function(expression, code, description) {
  result <- tryCatch(force(expression), error = identity)
  check(inherits(result, "ge_error") && identical(result$code, code), description)
}

work <- tempfile("pedigree-evaluation-tests-")
dir.create(work)
on.exit_cleanup <- function() {
  unlink(work, recursive = TRUE)
}

write_fixture <- function(data, name) {
  path <- file.path(work, name)
  write.csv(data, path, row.names = FALSE, na = "")
  path
}

tabular_source <- function(path, columns) {
  list(path = path, columns = columns, encoding = "UTF-8")
}

run_checks <- function() {
  check(identical(ge_id(c(" 00001 ", "1", "", NA_character_)),
                  c("00001", "1", NA_character_, NA_character_)),
        "Animal IDs retain leading zeros and normalize whitespace")

  check_close(ge_numeric(c("1.25D+01", "-999", ""), "trait", "-999"),
              c(12.5, NA_real_, NA_real_), "Numeric missing codes and D exponents")

  expect_error_code(ge_numeric("12 kg", "trait"), "ambiguous_numeric_values",
                    "Numeric columns reject embedded units")

  phenotype <- data.frame(
    animal = sprintf("%05d", seq_len(33L)),
    first = c(seq_len(30L), 300, 16, 81),
    second = c(seq_len(30L), 15, 400, 15),
    sex = rep(c("female", "male", "female"), 11L),
    stringsAsFactors = FALSE
  )
  phenotype$first[1L] <- NA_real_
  phenotype$second[2L] <- NA_real_
  phenotype <- rbind(phenotype,
                     data.frame(animal = "00034", first = NA_real_, second = NA_real_, sex = "female"),
                     data.frame(animal = "", first = 20, second = 20, sex = "male"))
  phenotype_path <- write_fixture(phenotype, "phenotypes.csv")

  analysis <- list(
    name = "synthetic_growth", species = "livestock", model = "animal",
    phenotypes = tabular_source(phenotype_path, list(id = "animal")),
    traits = list(list(name = "first", column = "first", unit = "unit"),
                  list(name = "second", column = "second", unit = "unit")),
    fixed = list(list(name = "sex", column = "sex", type = "class")),
    qc = list(method = "sd", sd = 3)
  )

  prepared <- ge_prepare_phenotypes(analysis, work)
  row31 <- match("00031", prepared$data$original_id)
  row32 <- match("00032", prepared$data$original_id)
  row33 <- match("00033", prepared$data$original_id)

  check(nrow(prepared$data) == 33L,
        "Partial-trait records remain after QC; missing IDs and empty records are omitted")
  check(is.na(prepared$data$y1[row31]) && is.finite(prepared$data$y2[row31]),
        "An outlier in the first trait does not erase the second trait")
  check(is.finite(prepared$data$y1[row32]) && is.na(prepared$data$y2[row32]),
        "An outlier in the second trait does not erase the first trait")
  check(prepared$data$y1[row33] == 81,
        "Three-SD QC runs once rather than recursively changing its limits")

  original_first <- phenotype$first[nzchar(phenotype$animal)]
  first_audit <- prepared$qc_audit[prepared$qc_audit$trait == "first", ]
  check_close(first_audit$upper, mean(original_first, na.rm = TRUE) +
                3 * sd(original_first, na.rm = TRUE),
              "QC limits come from the original eligible observations")
  check(identical(as.integer(prepared$qc_audit$removed), c(1L, 1L)),
        "Each trait reports its own rejected observation count")

  pedigree <- data.frame(
    animal = sprintf("%05d", seq_len(33L)),
    sire = "00090", dam = "00091", stringsAsFactors = FALSE
  )
  pedigree_path <- write_fixture(pedigree, "pedigree.csv")
  pedigree_source <- tabular_source(pedigree_path,
                                    list(id = "animal", sire = "sire", dam = "dam"))
  ped <- ge_prepare_pedigree(pedigree_source, prepared, character(), work)

  check(nrow(ped) == 35L && all(c("00090", "00091") %in% ped$original_id),
        "Referenced parents are included as identifiable founders")
  check(all(ped$sort[match(ped$sire[ped$sire > 0], ped$id)] <
              ped$sort[ped$sire > 0]), "Pedigree sorting places sires before offspring")
  check(all(ped$sort[match(ped$dam[ped$dam > 0], ped$id)] <
              ped$sort[ped$dam > 0]), "Pedigree sorting places dams before offspring")

  expect_error_code(ge_prepare_pedigree(pedigree_source, prepared, "not-in-pedigree", work),
                    "prediction_missing_pedigree", "Unresolved prediction targets are rejected")

  conflict <- rbind(pedigree, data.frame(animal = "00001", sire = "00092", dam = "00091"))
  conflict_source <- pedigree_source
  conflict_source$path <- write_fixture(conflict, "conflicting-pedigree.csv")
  expect_error_code(ge_prepare_pedigree(conflict_source, prepared, character(), work),
                    "conflicting_pedigree", "Contradictory parent records are rejected")

  loop <- data.frame(animal = c("00001", "00002"),
                     sire = c("00002", "00001"), dam = "0")
  loop_source <- pedigree_source
  loop_source$path <- write_fixture(loop, "cyclic-pedigree.csv")
  loop_phenotype <- list(data = data.frame(original_id = c("00001", "00002")))
  expect_error_code(ge_prepare_pedigree(loop_source, loop_phenotype, character(), work),
                    "pedigree_cycle", "Cycles in pedigree ancestry are rejected")

  counts <- data.frame(animal = c("A", "B", "C"), count = c(10, 11.5, 12))
  count_analysis <- analysis
  count_analysis$phenotypes <- tabular_source(write_fixture(counts, "counts.csv"),
                                             list(id = "animal"))
  count_analysis$traits <- list(list(name = "count", column = "count", unit = "count",
                                      distribution = "count_gaussian"))
  count_analysis$fixed <- list()
  expect_error_code(ge_prepare_phenotypes(count_analysis, work), "invalid_count_trait",
                    "Count traits reject fractional observations before estimation")

  repeated <- data.frame(animal = c("A", "A", "B", "B"),
                         date = c("2022-01-02", "2022/01/02", "2023-01-02", "2024-01-02"),
                         count = c(10, 11, 12, 13))
  repeated_analysis <- count_analysis
  repeated_analysis$model <- "repeatability"
  repeated_analysis$phenotypes$path <- write_fixture(repeated, "repeated.csv")
  repeated_analysis$event_date <- "date"
  repeated_analysis$record_key <- c("animal", "date")
  expect_error_code(ge_prepare_phenotypes(repeated_analysis, work), "duplicate_records",
                    "Equivalent date spellings cannot hide duplicate repeated records")

  numeric_date_source <- list(path = phenotype_path)
  expect_error_code(ge_date("45000", numeric_date_source, work), "ambiguous_date_origin",
                    "Numeric dates require an explicit or verified spreadsheet origin")
  numeric_date_source$date_origin <- "1899-12-30"
  check(identical(ge_date("45000", numeric_date_source, work),
                  as.Date(45000, origin = "1899-12-30")), "Declared date origin is used")

  traits <- data.frame(name = c("first", "second"), label = c("First", "Second"),
                       unit = c("unit", "unit"), stringsAsFactors = FALSE)
  components <- list(additive = matrix(c(4, -0.6, -0.6, 2), 2L),
                     litter = matrix(c(1, 0.2, 0.2, 0.5), 2L),
                     residual = matrix(c(5, 0.1, 0.1, 3), 2L))
  meta <- do.call(rbind, lapply(names(components), function(role) {
    data.frame(role = role, trait1 = c(1L, 1L, 2L), trait2 = c(1L, 2L, 2L))
  }))
  mixing <- diag(9L) * 0.05
  mixing[row(mixing) > col(mixing)] <- 0.008
  sampling_covariance <- tcrossprod(mixing)

  context <- list(nt = 2L, traits = traits, shared_role = "litter",
                  constraints = list(additive = matrix(FALSE, 2L, 2L),
                                     shared = matrix(FALSE, 2L, 2L),
                                     residual = matrix(FALSE, 2L, 2L)))
  fit <- list(components = list(additive = components$additive,
                                shared = components$litter,
                                residual = components$residual),
              param_meta = meta, param_cov = sampling_covariance, converged = TRUE)
  parameters <- ge_parameters(context, fit)

  check_close(parameters$estimate[parameters$parameter == "h2" & parameters$trait == "first"],
              0.4, "Heritability includes additive, litter and residual variance")

  parameter_values <- vapply(seq_len(nrow(meta)), function(index) {
    components[[meta$role[index]]][meta$trait1[index], meta$trait2[index]]
  }, numeric(1))

  evaluate <- function(values, kind, trait) {
    matrices <- components

    for (index in seq_along(values)) {
      role <- meta$role[index]
      first <- meta$trait1[index]
      second <- meta$trait2[index]
      matrices[[role]][first, second] <- values[index]
      matrices[[role]][second, first] <- values[index]
    }

    if (kind == "rg") {
      return(matrices$additive[1L, 2L] /
               sqrt(matrices$additive[1L, 1L] * matrices$additive[2L, 2L]))
    }

    total <- sum(vapply(matrices, function(value) value[trait, trait], numeric(1)))
    if (kind == "Vp") total else matrices$additive[trait, trait] / total
  }

  for (kind in c("h2", "Vp", "rg")) {
    for (trait in if (kind == "rg") 1L else 1:2) {
      numerical <- vapply(seq_along(parameter_values), function(index) {
        plus <- minus <- parameter_values
        plus[index] <- plus[index] + 1e-5
        minus[index] <- minus[index] - 1e-5
        (evaluate(plus, kind, trait) - evaluate(minus, kind, trait)) / 2e-5
      }, numeric(1))
      analytic <- ge_result_gradient(kind, components, meta, trait,
                                     if (kind == "rg") 2L else trait)
      check_close(analytic, numerical, paste(kind, "analytic gradient agrees with finite differences"),
                  tolerance = 1e-7)

      trait_name <- if (kind == "rg") "first : second" else traits$name[trait]
      estimate_se <- parameters$SE[parameters$parameter == kind & parameters$trait == trait_name]
      reference_se <- sqrt(drop(crossprod(numerical, sampling_covariance %*% numerical)))
      check_close(estimate_se, reference_se, paste(kind, "SE uses the full sampling covariance"),
                  tolerance = 1e-7)
    }
  }

  order <- c(9L, 2L, 5L, 1L, 6L, 4L, 8L, 3L, 7L)
  reordered <- fit
  reordered$param_meta <- meta[order, , drop = FALSE]
  reordered$param_cov <- sampling_covariance[order, order]
  check_close(ge_parameters(context, reordered)$SE, parameters$SE,
              "Changing backend parameter order does not change derived SE")

  missing_covariance <- fit
  missing_covariance$param_cov <- NULL
  missing_parameters <- ge_parameters(context, missing_covariance)
  check(all(is.na(missing_parameters$SE[missing_parameters$parameter %in% c("h2", "Vp", "rg")])),
        "Unavailable sampling covariance is reported instead of inventing SE")

  ebv <- data.frame(id = 1:3, trait = 1L, EBV = c(0.4, -0.2, 0.15),
                    SEP = c(1, 3, 1), F = c(0, 0, 0.25))
  reliability <- ge_result_reliability(context, fit, ebv)
  check_close(reliability$reliability[c(1L, 3L)], c(0.75, 0.8),
              "Reliability uses individual prior variance including inbreeding")
  check(is.na(reliability$reliability[2L]) &&
          grepl("not clipped", reliability$reliability_status[2L]),
        "An inconsistent PEV is flagged rather than clipped to zero")

  annual_context <- list(nt = 1L, traits = traits[1L, , drop = FALSE], trend_basis = "event_year",
                         data = data.frame(id = c(1L, 1L, 2L, 1L, 2L, 3L),
                                           event_year = c(2020L, 2020L, 2020L, 2021L, 2021L, 2022L),
                                           y1 = c(10, 12, 14, 13, 16, 16)))
  annual <- ge_annual(annual_context, ebv)$annual
  first_year <- annual[annual$year == 2020L, ]
  check(first_year$n_records == 3L && first_year$n_animals == 2L && first_year$n_ebv == 2L,
        "Repeated phenotypes do not multiply an animal's annual EBV contribution")
  check_close(first_year$mean_EBV, 0.1, "Annual mean EBV counts each animal once")
  check_close(first_year$mean_phenotype, 12, "Annual phenotype mean retains individual records")

  repeat_context <- list(nt = 1L, traits = traits[1L, , drop = FALSE], shared_role = "pe",
                        constraints = list(additive = matrix(FALSE, 1L, 1L),
                                           shared = matrix(FALSE, 1L, 1L),
                                           residual = matrix(FALSE, 1L, 1L)))
  repeat_fit <- list(components = list(additive = matrix(1, 1L, 1L),
                                       shared = matrix(2, 1L, 1L),
                                       residual = matrix(7, 1L, 1L)),
                    param_meta = data.frame(role = c("additive", "pe", "residual"),
                                            trait1 = 1L, trait2 = 1L),
                    param_cov = diag(c(0.1, 0.2, 0.3)), converged = TRUE)
  repeat_parameters <- ge_parameters(repeat_context, repeat_fit)
  check_close(repeat_parameters$estimate[repeat_parameters$parameter == "repeatability"],
              0.3, "Repeatability includes permanent environment in the numerator")
  check_close(repeat_parameters$estimate[repeat_parameters$parameter == "h2"],
              0.1, "Heritability excludes permanent environment from the numerator")

  ge_require(c("processx", "jsonlite"))
  script <- file.path(work, "subprocess check.R")
  writeLines(c("arguments <- commandArgs(trailingOnly = TRUE)",
               "stopifnot(identical(arguments, c('alpha beta', 'second')))",
               "cat(Sys.getenv('GE_TEST_ENV'), '\\n')"), script)
  extension <- if (.Platform$OS.type == "windows") ".exe" else ""
  executable <- file.path(R.home("bin"), paste0("Rscript", extension))
  process <- ge_run(executable, c("--vanilla", script, "alpha beta", "second"),
                    file.path(work, "subprocess"), env = c(GE_TEST_ENV = "environment-ok"))
  check(process$status == 0L &&
          identical(trimws(readLines(process$stdout)), "environment-ok"),
        "Process execution preserves spaced arguments and named environment values")

  cat(sprintf("Passed %d offline checks. No evaluation executable or private data was used.\n", checks))
}

tryCatch(run_checks(), finally = on.exit_cleanup())
