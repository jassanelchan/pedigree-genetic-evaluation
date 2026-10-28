run_general_data_checks <- function(check, check_close, expect_error_code, work, example_path) {
  specification <- ge_read_json(file.path(example_path, "multitrait-project.json"))
  analysis <- specification$analyses[[1L]]
  records <- ge_read_tabular(analysis$phenotypes, example_path)
  pedigree <- ge_read_tabular(specification$pedigree, example_path)
  fixture_number <- 0L

  write_records <- function(value, name = "records") {
    fixture_number <<- fixture_number + 1L
    path <- file.path(work, paste0("general-data-", fixture_number, "-", name, ".csv"))
    write.csv(value, path, row.names = FALSE, na = "")
    path
  }

  with_records <- function(value, definition = analysis) {
    definition$phenotypes$path <- write_records(value)
    definition
  }

  prepare <- function(definition = analysis, project = specification, base = example_path) {
    ge_prepare_context(project, definition, base, software = list())
  }

  fixed_term <- function(context, name) {
    context$fixed[[which(vapply(context$fixed, function(term) term$name == name, logical(1)))]]
  }

  context <- prepare()
  growth <- context$trait_rows[, "age120"]
  farrowing <- context$trait_rows[, "TNB"]
  sow_pe <- context$random$sow_pe
  litter <- context$random$birth_litter

  check(context$nt == 3L && nrow(context$data) == 1776L,
        "Mixed three-trait records are prepared without expanding growth measurements")
  check(nrow(pedigree) == 780L && nrow(context$ped) == 777L,
        "The pedigree retains the observed ancestral closure rather than unrelated founders")
  check(sum(growth) == 720L && sum(farrowing) == 1056L &&
          length(unique(context$data$id[farrowing])) == 352L,
        "Growth membership and repeated sow-event membership retain their distinct units")
  observed_counts <- colSums(!is.na(context$data[c("y1", "y2", "y3")]))
  check_close(observed_counts, c(720, 618, 1056),
              "Trait-specific observation counts survive mixed-record preparation")
  check(identical(dim(context$trait_rows), c(nrow(context$data), 3L)) &&
          all(context$trait_rows == as.matrix(context$data[paste0("member", 1:3)])),
        "Membership rows remain aligned with sorted data")
  check(all(context$trait_rows[, "bf120"] == growth) &&
          all(xor(growth, farrowing)) && any(growth & is.na(context$data$y2)),
        "Structural applicability is retained for a missing measurement")
  check(all(is.na(context$data$y1[farrowing])) && all(is.na(context$data$y2[farrowing])) &&
          all(is.na(context$data$y3[growth])),
        "A growth measurement is not copied onto the sow's later events")
  check(identical(names(context$random), c("additive", "birth_litter", "sow_pe")) &&
          identical(litter$traits, 1:2) && identical(sow_pe$traits, 3L),
        "Independent random blocks retain their explicit trait sets")
  check(all(context$data[[sow_pe$column]] == context$data$id) &&
          sow_pe$levels == nrow(context$ped), "Permanent environment is indexed by animal ID")
  check(length(unique(context$data$id[growth])) == sum(growth),
        "Each animal has at most one observed growth record")

  for (term_name in c("sex", "parity", "farrowing_season")) {
    term <- fixed_term(context, term_name)
    applicable <- rowSums(!is.na(context$data[paste0("y", term$traits)])) > 0L
    check(all(context$data[[term$column]][!applicable] == 1L),
          paste(term_name, "uses a neutral class code on inapplicable rows"))
  }

  covariate <- fixed_term(context, "measurement_covariate")
  check(all(context$data[[covariate$column]][is.na(context$data$y2)] == 0),
        "A trait-specific covariate uses zero only where that trait is unobserved")
  check(all(context$data[[litter$column]][farrowing] == 1L),
        "Missing birth-litter labels on farrowing-only rows receive a neutral placeholder")
  check(all(context$initial$birth_litter[3L, ] == 0) &&
          all(context$initial$sow_pe[-3L, ] == 0) &&
          all(context$constraints$birth_litter[3L, ]) &&
          all(context$constraints$sow_pe[-3L, ]),
        "Inactive random-component entries are structural zeros with matching constraints")
  check(isTRUE(context$constraints$residual[1L, 3L]) &&
          isTRUE(context$constraints$residual[2L, 3L]) &&
          context$initial$residual[1L, 3L] == 0 &&
          context$initial$residual[2L, 3L] == 0,
        "Unpaired residual covariances preserve the user's explicit zeros")

  for (term in context$random) {
    block <- context$initial[[term$name]][term$traits, term$traits, drop = FALSE]
    check(min(eigen(block, symmetric = TRUE, only.values = TRUE)$values) > 0,
          paste(term$name, "starts from a positive-definite active block"))
  }
  check(all(context$rank_checks$columns == context$rank_checks$rank),
        "Each trait checks only its applicable fixed-effect design")

  changed <- records
  no_backfat <- changed$event_type == "growth" & !nzchar(changed$bf120)
  changed$covariate[no_backfat] <- NA_character_
  missing_covariate <- prepare(with_records(changed))
  check(nrow(missing_covariate$data) == nrow(context$data),
        "A missing covariate is allowed when its only applicable trait was unmeasured")

  bad <- analysis
  bad$traits[[1L]]$record_filter <- NULL
  expect_error_code(prepare(bad), "missing_record_filter",
                    "Mixed single and repeated traits require explicit filters for every trait")

  changed <- records
  changed$age120[which(changed$event_type == "farrowing")[1L]] <- "185"
  expect_error_code(prepare(with_records(changed)), "observation_outside_record_filter",
                    "Values outside a trait's event filter are rejected rather than discarded")

  bad <- analysis
  bad$record_key <- NULL
  expect_error_code(prepare(bad), "missing_record_key",
                    "Repeated designs require an animal-and-event record key")
  bad$record_key <- "event_id"
  expect_error_code(prepare(bad), "missing_record_key",
                    "A repeated record key cannot omit the animal ID")
  bad$record_key <- "id"
  expect_error_code(prepare(bad), "missing_record_key",
                    "An animal ID alone does not distinguish repeated events")

  changed <- records
  changed$event_id[1L] <- NA_character_
  expect_error_code(prepare(with_records(changed)), "missing_record_key",
                    "An incomplete observed record key is rejected")
  changed <- rbind(records, records[1L, , drop = FALSE])
  expect_error_code(prepare(with_records(changed)), "duplicate_records",
                    "Duplicated animal-event records are rejected")

  bad <- analysis
  bad$traits[[1L]]$record_filter$values <- list("growth", "farrowing")
  changed <- records
  sow_id <- changed$id[changed$event_type == "farrowing"][1L]
  sow_growth <- changed$id == sow_id & changed$event_type == "growth"
  changed$age120[changed$id == sow_id] <- changed$age120[sow_growth]
  expect_error_code(prepare(with_records(changed, bad)), "duplicate_single_trait",
                    "Copying one growth phenotype onto multiple litter events is rejected")

  bad <- analysis
  bad$zero_covariances <- bad$zero_covariances[-1L]
  expect_error_code(prepare(bad), "unidentified_residual_covariance",
                    "An unpaired residual covariance needs an explicit zero constraint")

  bad <- analysis
  bad$random[[2L]]$column <- "id"
  expect_error_code(prepare(bad), "invalid_pe_group",
                    "Permanent environment does not accept a custom grouping column")

  farrow_rows <- records$event_type == "farrowing"
  retain <- !farrow_rows
  retain[which(farrow_rows)[!duplicated(records$id[farrow_rows])]] <- TRUE
  expect_error_code(prepare(with_records(records[retain, , drop = FALSE])), "unidentified_pe",
                    "PE is not estimated when every sow has only one observed event")

  changed <- records
  changed$unique_group <- seq_len(nrow(changed))
  bad <- analysis
  bad$random <- c(bad$random, list(list(name = "event_group", role = "iid",
                                      column = "unique_group", traits = list("TNB"))))
  expect_error_code(prepare(with_records(changed, bad)), "unidentified_random_effect",
                    "An independent group seen once per observation is confounded with residual")

  bad <- analysis
  bad$random <- c(bad$random, list(list(name = "same_litter", role = "iid",
                                      column = "litter", traits = list("age120"))))
  expect_error_code(prepare(bad), "confounded_random_effects",
                    "Overlapping trait sets cannot hide identical random-effect partitions")

  bad <- analysis
  bad$shared <- list(column = "litter")
  expect_error_code(prepare(bad), "inapplicable_shared_definition",
                    "The general model cannot silently ignore a legacy shared definition")

  bad <- analysis
  bad$inbreeding <- "true"
  expect_error_code(prepare(bad), "invalid_boolean",
                    "A string cannot silently disable pedigree inbreeding")
  bad <- analysis
  bad$allow_observed_founders <- "false"
  expect_error_code(prepare(bad), "invalid_boolean",
                    "The observed-founder policy requires a real Boolean")
  bad_source <- analysis$phenotypes
  bad_source$header <- "true"
  expect_error_code(ge_read_tabular(bad_source, example_path), "invalid_boolean",
                    "The table-header option cannot be a string Boolean")

  for (enabled in c(TRUE, FALSE)) {
    chosen <- analysis
    chosen$inbreeding <- enabled
    chosen$allow_observed_founders <- enabled
    chosen$phenotypes$header <- TRUE
    chosen_context <- prepare(chosen)
    check(identical(chosen_context$inbreeding, enabled),
          paste("An explicit logical inbreeding option remains", enabled))
  }
  headerless_source <- analysis$phenotypes
  headerless_source$header <- FALSE
  headerless <- ge_read_tabular(headerless_source, example_path)
  check(nrow(headerless) == nrow(records) + 1L && headerless[[1L]][1L] == "id",
        "A real FALSE header option preserves the first row as data")

  bad <- analysis
  bad$random[[1L]]$traits <- list("not_a_trait")
  expect_error_code(prepare(bad), "invalid_effect_traits",
                    "Random effects require known trait names")
  bad <- analysis
  bad$fixed[[1L]]$traits <- list("age120", "age120")
  expect_error_code(prepare(bad), "invalid_effect_traits",
                    "Fixed-effect active trait names must be distinct")

  changed <- records
  changed$litter[which(changed$event_type == "growth")[1L]] <- NA_character_
  expect_error_code(prepare(with_records(changed)), "missing_litter",
                    "An applicable litter grouping cannot be missing")
  changed <- records
  changed$sex[which(changed$event_type == "growth")[1L]] <- NA_character_
  expect_error_code(prepare(with_records(changed)), "missing_fixed_effect",
                    "An applicable observed fixed effect cannot be missing")

  bad <- analysis
  bad$initial <- list(additive = diag(3L))
  expect_error_code(prepare(bad), "invalid_initial_matrix",
                    "A zero free covariance is not accepted as a starting value")
  bad$initial$additive <- matrix(c(2, 0.1, 0.1, 0.2, 2, 0.1, 0.1, 0.1, 2), 3L, byrow = TRUE)
  expect_error_code(prepare(bad), "invalid_initial_matrix",
                    "Asymmetric starting matrices are rejected")
  bad$initial$additive <- matrix(2, 3L, 3L)
  diag(bad$initial$additive) <- 1
  expect_error_code(prepare(bad), "invalid_initial_matrix",
                    "A non-positive-definite starting matrix is rejected")
  bad$initial <- list(birth_litter = matrix(c(1, 0.1, 0, 0.1, 1, 0, 0, 0, 1), 3L))
  expect_error_code(prepare(bad), "invalid_initial_matrix",
                    "Starting values cannot populate an inactive trait block")
  bad$initial <- list(additive = list(list(1, 0.1), list(0.1, 1)))
  expect_error_code(prepare(bad), "invalid_initial_matrix",
                    "Starting matrices must use the complete global trait dimensions")

  for (index in list(1.5, -1L, 0L, 4L, NA_real_)) {
    bad <- analysis
    bad$zero_covariances[[1L]]$trait1 <- index
    expect_error_code(prepare(bad), "unsupported_constraint",
                      paste("Zero constraints reject an invalid trait index:", index))
  }
  bad <- analysis
  bad$zero_covariances[[1L]]$trait2 <- 1L
  expect_error_code(prepare(bad), "unsupported_constraint",
                    "Zero-covariance constraints cannot remove an active variance")
  bad <- analysis
  bad$zero_covariances[[1L]]$component <- "unknown_component"
  expect_error_code(prepare(bad), "unsupported_constraint",
                    "Zero constraints require an existing component")

  changed <- records
  changed$disjoint_group <- ifelse(changed$event_type == "growth",
                                    paste0("growth_", changed$litter),
                                    paste0("sow_", changed$id))
  bad <- analysis
  bad$random <- list(list(name = "event_block", role = "iid", column = "disjoint_group",
                          traits = list("age120", "TNB")))
  expect_error_code(prepare(with_records(changed, bad)), "unidentified_random_covariance",
                    "An identity-term covariance needs groups represented in both active traits")
  bad$zero_covariances <- c(bad$zero_covariances,
                            list(list(component = "event_block", trait1 = 1L, trait2 = 3L)))
  constrained <- prepare(with_records(changed, bad))
  check(isTRUE(constrained$constraints$event_block[1L, 3L]) &&
          constrained$initial$event_block[1L, 3L] == 0,
        "Explicitly fixing a disconnected identity covariance preserves the otherwise valid model")

  many_records <- data.frame(id = sprintf("G%04d", seq_len(30L)), stringsAsFactors = FALSE)
  for (trait in seq_len(5L)) {
    many_records[[paste0("trait", trait)]] <- 10 * trait +
      sin(seq_len(nrow(many_records)) * trait) + seq_len(nrow(many_records)) / (trait + 1)
  }
  many_pedigree <- data.frame(id = many_records$id, sire = "parent1", dam = "parent2")
  many_spec <- list(pedigree = list(path = write_records(many_pedigree, "five-trait-pedigree"),
                                   columns = list(id = "id", sire = "sire", dam = "dam")))
  many_analysis <- list(name = "five_traits", species = "livestock", model = "animal",
                        phenotypes = list(path = write_records(many_records, "five-trait-records"),
                                          columns = list(id = "id")),
                        traits = lapply(seq_len(5L), function(trait) {
                          list(name = paste0("trait", trait), column = paste0("trait", trait), unit = "unit")
                        }),
                        qc = list(method = "none"))
  for (family in c("animal", "general")) {
    many_analysis$model <- family
    many_context <- prepare(many_analysis, many_spec)
    check(many_context$nt == 5L && nrow(many_context$data) == 30L &&
            identical(dim(many_context$initial$additive), c(5L, 5L)),
          paste("Five traits have no artificial two- or three-trait cap in", family, "models"))
  }

  repeated_records <- data.frame(id = rep(sprintf("R%03d", 1:8), each = 3L),
                                  event = rep(1:3, 8L),
                                  first = seq_len(24L) + sin(seq_len(24L)),
                                  second = 10 + cos(seq_len(24L)) + seq_len(24L) / 4)
  repeated_pedigree <- data.frame(id = unique(repeated_records$id), sire = "father", dam = "mother")
  repeated_spec <- list(pedigree = list(path = write_records(repeated_pedigree, "repeated-pedigree"),
                                       columns = list(id = "id", sire = "sire", dam = "dam")))
  repeated_analysis <- list(name = "two_repeated_traits", species = "livestock", model = "repeatability",
                            phenotypes = list(path = write_records(repeated_records, "repeated-records"),
                                              columns = list(id = "id")),
                            traits = list(list(name = "first", column = "first", unit = "unit"),
                                          list(name = "second", column = "second", unit = "unit")),
                            record_key = c("id", "event"), qc = list(method = "none"))
  repeated_context <- prepare(repeated_analysis, repeated_spec)
  check(repeated_context$nt == 2L && nrow(repeated_context$data) == 24L &&
          all(repeated_context$traits$record_type == "repeated") &&
          identical(repeated_context$random$shared$traits, 1:2),
        "The repeatability shortcut supports two repeated traits on one PE block")
  check(all(repeated_context$data$shared == repeated_context$data$id) &&
          repeated_context$random$shared$levels == nrow(repeated_context$ped),
        "Multitrait PE remains animal-indexed after pedigree recoding")

  invisible(context)
}
