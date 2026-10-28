args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
example_dir <- if (length(args)) args[[1L]] else dirname(
  normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/")
)
dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite first.")
set.seed(20261028)
draw <- function(n, covariance) {
  matrix(rnorm(n * nrow(covariance)), n, nrow(covariance)) %*% chol(covariance)
}
n_founders <- 60L
n_litters <- 120L
n_animals <- n_litters * 6L
sire <- rep(sample(seq_len(20L), n_litters, TRUE), each = 6L)
dam <- rep(sample(21:60, n_litters, TRUE), each = 6L)
ids <- sprintf("SYN_%04d", seq_len(n_founders + n_animals))
additive_covariance <- matrix(c(4, 0.6, 0.4, 0.6, 2, 0.2, 0.4, 0.2, 1.5), 3L)
litter_covariance <- matrix(c(1, 0.1, 0.1, 0.5), 2L)
residual_covariance <- matrix(c(6, 0.5, 0.5, 3), 2L)
founders <- draw(n_founders, additive_covariance)
additive <- 0.5 * founders[sire, ] + 0.5 * founders[dam, ] +
  draw(n_animals, 0.5 * additive_covariance)
litter <- draw(n_litters, litter_covariance)[rep(seq_len(n_litters), each = 6L), ]
permanent <- rnorm(n_animals, sd = sqrt(1.2))
sex <- sample(c("female", "male"), n_animals, TRUE)
covariate <- rnorm(n_animals)
growth <- additive[, 1:2] + litter + draw(n_animals, residual_covariance) +
  cbind(180 + 2 * (sex == "male"), 12 + 0.4 * covariate)
growth[seq(7L, n_animals, by = 7L), 2L] <- NA_real_
sows <- which(sex == "female")
who <- rep(sows, each = 3L)
parity <- rep(seq_len(3L), length(sows))
season <- sample(seq_len(4L), length(who), TRUE)
total_born <- round(15 + additive[who, 3L] + permanent[who] +
                      rnorm(length(who), sd = 2) + 0.5 * parity + 0.2 * season)
stopifnot(all(total_born >= 0))
birth <- as.Date("2017-01-01") + (seq_len(n_animals) %% 6L) * 365L
growth_records <- data.frame(
  id = tail(ids, n_animals), event_id = "growth", event_type = "growth",
  age120 = growth[, 1L], bf120 = growth[, 2L], TNB = NA_real_,
  sex = sex, covariate = covariate, parity = NA_integer_, season = NA_integer_,
  litter = sprintf("L_%03d", rep(seq_len(n_litters), each = 6L)),
  event_date = as.character(birth + 180L), stringsAsFactors = FALSE)
repeated_records <- data.frame(
  id = tail(ids, n_animals)[who], event_id = paste0("farrowing_", parity),
  event_type = "farrowing", age120 = NA_real_, bf120 = NA_real_, TNB = total_born,
  sex = NA_character_, covariate = NA_real_, parity = parity, season = season,
  litter = NA_character_, event_date = as.character(birth[who] + 400L + 190L * (parity - 1L)),
  stringsAsFactors = FALSE)
records <- rbind(growth_records, repeated_records)
pedigree <- data.frame(
  id = ids, sire = c(rep("0", n_founders), ids[sire]),
  dam = c(rep("0", n_founders), ids[dam]), stringsAsFactors = FALSE)
write.csv(records, file.path(example_dir, "mixed-records.csv"), row.names = FALSE, na = "")
write.csv(pedigree, file.path(example_dir, "mixed-pedigree.csv"), row.names = FALSE, na = "")
trait <- function(name, unit, record_type, event, distribution = "gaussian") {
  list(name = name, column = name, unit = unit, distribution = distribution,
       record_type = record_type, record_filter = list(column = "event_type", values = list(event)))
}
project <- list(
  schema_version = 2L, engine = "compare", output_dir = "example_runs",
  timeout = 300L, threads = 2L,
  pedigree = list(path = "mixed-pedigree.csv",
                  columns = list(id = "id", sire = "sire", dam = "dam")),
  analyses = list(list(
    name = "mixed_growth_and_litter_size", species = "pig", model = "general",
    phenotypes = list(path = "mixed-records.csv", columns = list(id = "id")),
    traits = list(trait("age120", "days", "single", "growth"),
                  trait("bf120", "mm", "single", "growth"),
                  trait("TNB", "piglets/litter", "repeated", "farrowing", "count_gaussian")),
    fixed = list(
      list(name = "sex", column = "sex", traits = list("age120", "bf120")),
      list(name = "measurement_covariate", column = "covariate", type = "cov",
           traits = list("bf120")),
      list(name = "parity", column = "parity", traits = list("TNB")),
      list(name = "farrowing_season", column = "season", traits = list("TNB"))),
    random = list(
      list(name = "birth_litter", role = "litter", column = "litter",
           traits = list("age120", "bf120")),
      list(name = "sow_pe", role = "pe", traits = list("TNB"))),
    record_key = list("id", "event_id"), event_date = "event_date",
    qc = list(method = "none"),
    zero_covariances = list(
      list(component = "residual", trait1 = 1L, trait2 = 3L),
      list(component = "residual", trait1 = 2L, trait2 = 3L)))))
jsonlite::write_json(project, file.path(example_dir, "multitrait-project.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Synthetic example:", nrow(records), "records,", nrow(pedigree),
    "pedigree animals,", length(sows), "recorded sows.\n")
