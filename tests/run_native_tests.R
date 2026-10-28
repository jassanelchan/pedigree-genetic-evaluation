arguments <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  at <- match(name, arguments)
  if (is.na(at)) return(default)
  if (at == length(arguments)) stop("Missing value after ", name)
  arguments[[at + 1L]]
}
profile <- option("--profile")
if (is.null(profile)) {
  stop("Usage: Rscript tests/run_native_tests.R --profile host_profile.json [--output DIRECTORY]")
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
repository <- dirname(dirname(normalizePath(sub("^--file=", "", script_argument[[1L]]),
                                             winslash = "/", mustWork = TRUE)))
skill <- file.path(repository, "skills", "pedigree-genetic-evaluation")
source(file.path(skill, "scripts", "common.R"), encoding = "UTF-8")
ge_require(c("jsonlite", "digest", "processx"))
software <- ge_read_json(ge_path(profile))$software
if (!length(software)) stop("The profile must contain software executable paths.")
output <- option("--output", tempfile("pedigree-native-"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- ge_path(output)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
example <- file.path(skill, "references", "example")
entry <- file.path(skill, "scripts", "run_eval.R")
reports <- list()

run_case <- function(filename, label) {
  case_directory <- ge_stage_directory(output, label)
  project <- ge_read_json(file.path(example, filename))
  project$engine <- "compare"
  project$software <- software
  project$output_dir <- "runs"
  project$pedigree$path <- ge_path(project$pedigree$path, example)
  for (index in seq_along(project$analyses)) {
    project$analyses[[index]]$phenotypes$path <-
      ge_path(project$analyses[[index]]$phenotypes$path, example)
  }
  project_file <- file.path(case_directory, "project.json")
  ge_write_json(project, project_file)
  process <- ge_run(rscript, c("--vanilla", entry, "--spec", project_file),
                    case_directory, log_name = "evaluation.stdout.log", timeout = 2400)
  if (process$status != 0L) stop("Native run failed; inspect ", process$stdout)
  runs <- list.dirs(file.path(case_directory, "runs"), recursive = FALSE, full.names = TRUE)
  if (length(runs) != 1L) stop("Expected one fresh run directory.")
  run <- runs[[1L]]
  manifest <- ge_read_json(file.path(run, "run_manifest.json"))
  if (!identical(manifest$state, "complete") || isTRUE(manifest$comparison_failed)) {
    stop("Incomplete evaluation or backend discrepancy: ", run)
  }
  comparisons <- list.files(file.path(run, "delivery"), pattern = "^comparison\\.json$",
                             recursive = TRUE, full.names = TRUE)
  if (length(comparisons) != length(project$analyses)) stop("Missing comparison reports.")
  for (path in comparisons) {
    compared <- ge_read_json(path)
    if (!isTRUE(compared$pass) || !isTRUE(compared$fitted_verified)) {
      stop("Component, EBV/PEV or fitted-value comparison failed: ", path)
    }
  }
  logs <- list.files(run, pattern = "\\.process\\.json$", recursive = TRUE, full.names = TRUE)
  before <- setNames(vapply(logs, ge_hash_file, character(1)), logs)
  resumed <- ge_run(rscript, c("--vanilla", entry, "--spec", project_file, "--resume", run),
                    case_directory, log_name = "resume.stdout.log", timeout = 2400)
  after_logs <- list.files(run, pattern = "\\.process\\.json$", recursive = TRUE, full.names = TRUE)
  after <- setNames(vapply(after_logs, ge_hash_file, character(1)), after_logs)
  if (resumed$status != 0L || !identical(before, after)) {
    stop("Resume did not reuse the verified native results: ", run)
  }
  text <- readLines(resumed$stdout, warn = FALSE)
  if (!any(grepl("Reusing verified fit:", text, fixed = TRUE)) ||
      !any(grepl("Reusing verified comparison:", text, fixed = TRUE))) {
    stop("Expected explicit fit and comparison cache reuse.")
  }
  list(state = "pass", project = filename, run_directory = run,
       native_processes = length(logs), comparisons = length(comparisons),
       resume_reused = TRUE)
}
reports$legacy <- run_case("project.json", "legacy")
reports$general <- run_case("multitrait-project.json", "general")
ge_write_json(reports, file.path(output, "native_test_report.json"))
cat("Native BLUPF90/DMU evaluation, comparison and resume checks passed.\n")
cat("Report:", file.path(output, "native_test_report.json"), "\n")
