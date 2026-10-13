ge_stop <- function(message, code = "invalid_input") {
  stop(structure(list(message = message, call = NULL, code = code),
                 class = c("ge_error", "error", "condition")))
}

ge_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) ge_stop(paste("Missing R packages:", paste(missing, collapse = ", ")),
                               "missing_dependency")
}

ge_write_json <- function(value, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, na = "null", null = "null")
  invisible(file)
}

ge_read_json <- function(file) {
  jsonlite::fromJSON(file, simplifyVector = FALSE)
}

ge_matrix_lines <- function(x) {
  apply(as.matrix(x), 1L, function(row) paste(sprintf("%.17g", row), collapse = " "))
}

ge_hash_file <- function(file) {
  if (!file.exists(file)) ge_stop(paste("File not found:", file), "missing_file")
  digest::digest(file = file, algo = "sha256")
}

ge_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

ge_path <- function(x, base = getwd(), must_work = TRUE) {
  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    ge_stop("A nonempty file path is required", "missing_path")
  }
  absolute <- grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", x)
  if (!absolute) x <- file.path(base, x)
  normalizePath(x, winslash = "/", mustWork = must_work)
}

ge_id <- function(x) {
  value <- trimws(as.character(x))
  value[is.na(x) | toupper(value) %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  value
}

ge_parent_id <- function(x) {
  value <- ge_id(x)
  value[!is.na(value) & value %in% c("0", "-", ".")] <- NA_character_
  value
}

ge_numeric <- function(x, field, missing = character()) {
  value <- trimws(as.character(x))
  absent <- is.na(x) | toupper(value) %in% c("", "NA", "N/A", "NULL", "NAN") |
    value %in% as.character(missing)
  parsed <- suppressWarnings(as.numeric(gsub("[dD]", "E", value)))
  invalid <- !absent & (!is.finite(parsed) | is.na(parsed))
  if (any(invalid)) {
    examples <- unique(value[invalid])[seq_len(min(5L, length(unique(value[invalid]))))]
    ge_stop(paste0("Numeric conversion failed in ", field, ": ", paste(examples, collapse = ", ")),
            "ambiguous_numeric_values")
  }
  parsed[absent] <- NA_real_
  parsed
}

ge_column <- function(data, column, role) {
  if (is.null(column)) ge_stop(paste("Column required for", role), "missing_column_mapping")
  if (is.numeric(column) || grepl("^[0-9]+$", as.character(column))) {
    index <- as.integer(column)
    if (index < 1L || index > ncol(data)) ge_stop(paste("Column out of range:", role))
    return(data[[index]])
  }
  if (!column %in% names(data)) ge_stop(paste("Column not found:", column, "for", role),
                                      "missing_column")
  data[[column]]
}

ge_read_tabular <- function(source, base = getwd()) {
  path <- ge_path(source$path, base)
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("xlsx", "xls")) {
    ge_require("readxl")
    sheets <- readxl::excel_sheets(path)
    sheet <- source$sheet
    if (is.null(sheet)) {
      if (length(sheets) != 1L) ge_stop(paste("Choose a sheet for", path,
                                            "Available:", paste(sheets, collapse = ", ")),
                                      "ambiguous_sheet")
      sheet <- sheets[[1L]]
    }
    data <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text",
                                           .name_repair = "minimal"), check.names = FALSE)
  } else {
    sep <- source$sep
    if (is.null(sep)) sep <- if (extension == "csv") "," else if (extension == "tsv") "\t" else ""
    header <- if (is.null(source$header)) TRUE else isTRUE(source$header)
    encoding <- if (is.null(source$encoding)) "UTF-8-BOM" else source$encoding
    data <- read.table(path, sep = sep, header = header, colClasses = "character",
                       check.names = FALSE, stringsAsFactors = FALSE, comment.char = "",
                       quote = if (sep == "") "" else "\"", fileEncoding = encoding,
                       na.strings = character(), fill = FALSE)
  }
  if (anyDuplicated(names(data))) ge_stop(paste("Duplicate column names in", path),
                                         "ambiguous_column_names")
  attr(data, "source_path") <- path
  data
}

ge_date <- function(x, source, base = getwd()) {
  value <- trimws(as.character(x))
  absent <- is.na(x) | toupper(value) %in% c("", "NA", "N/A", "NULL")
  numeric_rows <- !absent & grepl("^[0-9]+(\\.[0-9]+)?$", value)
  result <- as.Date(rep(NA_character_, length(value)))
  if (any(numeric_rows)) {
    origin <- source$date_origin
    path <- ge_path(source$path, base)
    if (is.null(origin) && tolower(tools::file_ext(path)) == "xlsx") {
      ge_require("openxlsx")
      origin <- openxlsx::getDateOrigin(path)
      if (origin == "1900-01-01") origin <- "1899-12-30"
    }
    if (is.null(origin)) ge_stop("Numeric dates need a verified Excel date origin",
                                "ambiguous_date_origin")
    result[numeric_rows] <- as.Date(as.numeric(value[numeric_rows]), origin = origin)
  }
  text_rows <- which(!absent & !numeric_rows)
  if (length(text_rows)) {
    format <- source$date_format
    if (is.null(format)) {
      cleaned <- gsub("/", "-", substr(value[text_rows], 1L, 10L))
      if (any(!grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", cleaned))) {
        ge_stop("Text dates require an unambiguous YYYY-MM-DD format or explicit date_format",
                "ambiguous_date_format")
      }
      result[text_rows] <- as.Date(cleaned, format = "%Y-%m-%d")
    } else result[text_rows] <- as.Date(value[text_rows], format = format)
    if (anyNA(result[text_rows])) ge_stop("Invalid date values", "invalid_date")
  }
  result
}

ge_run <- function(executable, args = character(), directory, stdin = NULL,
                   log_name = "stdout.log", timeout = 7200, env = character()) {
  executable <- ge_path(executable)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- ge_path(directory)
  ge_require("processx")
  stdout <- file.path(directory, log_name)
  stderr <- paste0(stdout, ".stderr")
  input <- ""
  if (!is.null(stdin)) {
    input <- file.path(directory, "stdin.txt")
    writeLines(as.character(stdin), input, useBytes = TRUE)
  }
  command <- paste(shQuote(executable), paste(shQuote(args), collapse = " "))
  environment <- Sys.getenv()
  if (length(env)) {
    if (is.null(names(env)) || any(!nzchar(names(env)))) {
      pairs <- strsplit(env, "=", fixed = TRUE)
      if (any(lengths(pairs) < 2L)) ge_stop("Process env entries must be NAME=value")
      names(env) <- vapply(pairs, `[[`, character(1), 1L)
      env <- setNames(vapply(pairs, function(x) paste(x[-1L], collapse = "="), character(1)), names(env))
    }
    environment[names(env)] <- env
  }
  start <- Sys.time()
  process <- tryCatch(processx::run(executable, args = args, wd = directory,
                                    stdout = stdout, stderr = stderr,
                                    stdin = if (nzchar(input)) input else NULL,
                                    timeout = timeout, env = environment,
                                    error_on_status = FALSE, cleanup_tree = TRUE,
                                    windows_hide_window = TRUE), error = function(e) e)
  status <- if (inherits(process, "error")) {
    message <- conditionMessage(process)
    cat(message, "\n", file = stderr, append = TRUE)
    if (grepl("timed out|timeout", message, ignore.case = TRUE)) 124L else 1L
  } else if (isTRUE(process$timeout)) {
    cat("Process timed out; its process tree was terminated.\n", file = stderr, append = TRUE)
    124L
  } else process$status
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  result <- list(status = as.integer(status), stdout = stdout, stderr = stderr,
                 elapsed = elapsed, command = command)
  ge_write_json(result, paste0(stdout, ".process.json"))
  result
}

ge_stage_directory <- function(root, name) {
  directory <- file.path(root, name)
  if (dir.exists(directory)) {
    number <- 2L
    while (dir.exists(paste0(directory, "_", number))) number <- number + 1L
    directory <- paste0(directory, "_", number)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ge_path(directory)
}

ge_write_bom_csv <- function(x, path) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.raw(c(239, 187, 191)), connection)
  lines <- capture.output(write.csv(x, row.names = FALSE, na = ""))
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\r\n"), "\r\n"))), connection)
}
