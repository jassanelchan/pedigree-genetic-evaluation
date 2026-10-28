run_general_results_checks <- function(check, check_close, expect_error_code, work) {
  output <- file.path(work, "general-results")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)

  term <- function(name, role, traits) list(name = name, role = role, traits = as.integer(traits))
  random <- list(term("additive", "additive", 1:4), term("birth_litter", "litter", 1:2),
                 term("nurse_litter", "litter", 1:2), term("pe_tnb", "pe", 3L),
                 term("pe_lactation", "pe", 4L), term("pen", "iid", c(1, 3, 4)))

  mat <- function(active, diagonal, cross = 0.05) {
    value <- matrix(0, 4, 4)
    value[active, active] <- cross
    diag(value)[active] <- diagonal
    value
  }

  components <- list(additive = mat(1:4, c(4, 2, 1, 0.8), 0.15),
                     birth_litter = mat(1:2, c(1, 0.5)),
                     nurse_litter = mat(1:2, c(0.7, 0.3), -0.04),
                     pe_tnb = mat(3L, 0.6), pe_lactation = mat(4L, 0.9),
                     pen = mat(c(1, 3, 4), c(0.4, 0.3, 0.2)),
                     residual = mat(1:4, c(5, 3, 7, 2), 0.2))
  components$additive[1, 4] <- components$additive[4, 1] <- 0
  components$residual[1:2, 3:4] <- components$residual[3:4, 1:2] <- 0

  ctx <- list(nt = 4L, family = "general", name = "anonymous_general", species = "livestock",
              traits = data.frame(name = c("age", "bf", "TNB", "milk"),
                                   label = c("Age", "Backfat", "Litter size", "Milk"),
                                   unit = c("days", "mm", "piglets/litter", "kg")),
              random = random, trend_basis = "event_year", shared_role = NULL,
              model_description = "Four-trait animal model with two litter terms, two disjoint PE terms and pen.")
  ctx$constraints <- lapply(components, function(x) {
    key <- names(components)[which(vapply(components, identical, logical(1), x))[1]]
    index <- match(key, vapply(random, function(x) x$name, character(1)))
    active <- if (key == "residual") 1:4 else random[[index]]$traits
    answer <- matrix(TRUE, 4, 4)
    answer[active, active] <- FALSE
    answer
  })
  ctx$constraints$additive[1, 4] <- ctx$constraints$additive[4, 1] <- TRUE
  ctx$constraints$residual[1:2, 3:4] <- ctx$constraints$residual[3:4, 1:2] <- TRUE

  meta <- do.call(rbind, lapply(names(components), function(key) {
    indices <- which(upper.tri(ctx$constraints[[key]], diag = TRUE) & !ctx$constraints[[key]], arr.ind = TRUE)
    data.frame(role = key, component = key, trait1 = indices[, 1], trait2 = indices[, 2])
  }))
  mix <- diag(seq(0.05, 0.10, length.out = nrow(meta)))
  mix[lower.tri(mix)] <- 0.003
  covariance <- tcrossprod(mix)
  fit <- list(components = components, param_meta = meta, param_cov = covariance,
              converged = TRUE, engine = "anonymous", version = "test", rounds = 1L,
              minus2logL = 100, warnings = character())
  parameters <- ge_parameters(ctx, fit)
  canonical <- ge_result_components(ctx, fit)

  evaluate <- function(values, kind, first, second = first) {
    matrices <- components
    for (index in seq_len(nrow(meta))) {
      key <- meta$role[index]
      one <- meta$trait1[index]
      two <- meta$trait2[index]
      matrices[[key]][one, two] <- matrices[[key]][two, one] <- values[index]
    }

    total <- sum(vapply(matrices, function(x) x[first, first], numeric(1)))
    additive <- matrices$additive[first, first]
    persistent <- names(random)[0]
    persistent <- vapply(Filter(function(x) x$role == "pe" && first %in% x$traits, random),
                          function(x) x$name, character(1))

    switch(kind, Vp = total, h2 = additive / total,
           repeatability = (additive + sum(vapply(matrices[persistent], function(x) x[first, first], numeric(1)))) / total,
           rg = matrices$additive[first, second] / sqrt(additive * matrices$additive[second, second]))
  }

  values <- vapply(seq_len(nrow(meta)), function(i) components[[meta$role[i]]][meta$trait1[i], meta$trait2[i]], numeric(1))

  for (kind in c("Vp", "h2", "repeatability", "rg")) {
    pairs <- if (kind == "rg") t(combn(1:4, 2)) else cbind(1:4, 1:4)
    if (kind == "repeatability") pairs <- pairs[3:4, , drop = FALSE]

    for (index in seq_len(nrow(pairs))) {
      first <- pairs[index, 1]
      second <- pairs[index, 2]
      numerical <- vapply(seq_along(values), function(i) {
        step <- 1e-5 * max(1, abs(values[i]))
        plus <- minus <- values
        plus[i] <- plus[i] + step
        minus[i] <- minus[i] - step
        (evaluate(plus, kind, first, second) - evaluate(minus, kind, first, second)) / (2 * step)
      }, numeric(1))
      analytic <- ge_result_gradient(kind, canonical, meta, first, second)
      check(max(abs(numerical - analytic)) < 1e-8, paste("Derivative", kind, first, second))

      estimate <- parameters[parameters$parameter == kind & parameters$trait1 == first & parameters$trait2 == second, ]
      se <- sqrt(drop(crossprod(numerical, covariance %*% numerical)))
      check(nrow(estimate) == 1L && abs(estimate$SE - se) < 1e-8, paste("Full delta SE", kind, first, second))
    }
  }

  check(sum(parameters$parameter == "rg") == 6L, "All six genetic correlations present")
  check(sum(parameters$parameter == "repeatability") == 2L, "Only PE-applicable traits have repeatability")
  check(abs(subset(parameters, parameter == "Vp" & trait == "age")$estimate - 11.1) < 1e-12,
        "Both litter terms and iid pen enter age variance")
  check(abs(subset(parameters, parameter == "repeatability" & trait == "TNB")$estimate - 1.6 / 8.9) < 1e-12,
        "PE enters repeatability; iid does not")
  check(all(parameters$estimate[parameters$SE_method == "inactive"] == 0) &&
          all(is.na(parameters$SE[parameters$SE_method == "inactive"])), "Inactive structural zeros have no fake sampling SE")
  check(subset(parameters, parameter == "additive_covariance" & trait1 == 1 & trait2 == 4)$SE_method == "fixed",
        "Explicit zero separated from inactive zero")
  check(subset(parameters, parameter == "rg" & trait1 == 1 & trait2 == 4)$SE == 0,
        "Constrained-zero genetic correlation has zero Delta SE")

  set.seed(42)
  order <- sample(nrow(meta))
  reordered <- fit
  reordered$param_meta <- meta[order, , drop = FALSE]
  reordered$param_cov <- covariance[order, order]
  check(isTRUE(all.equal(ge_parameters(ctx, reordered), parameters)), "Arbitrary sampling order is invariant")
  incomplete <- fit
  incomplete$param_meta <- meta[-1, ]
  incomplete$param_cov <- covariance[-1, -1]
  check(all(is.na(subset(ge_parameters(ctx, incomplete), parameter %in% c("h2", "Vp", "rg"))$SE)),
        "Incomplete sampling metadata never generates partial SE")
  unknown <- fit
  unknown$components$pe_tnb[1, 1] <- 1
  check(inherits(try(ge_parameters(ctx, unknown), silent = TRUE), "try-error"), "Inactive nonzero component rejected")

  fractional <- fit
  fractional$param_meta$trait1[1] <- 1.5
  check(all(is.na(subset(ge_parameters(ctx, fractional), parameter == "h2")$SE)),
        "Noninteger sampling trait indices are rejected rather than truncated")
  extra <- fit
  extra$components$unmodeled <- diag(4)
  check(inherits(try(ge_parameters(ctx, extra), silent = TRUE), "try-error"),
        "Undeclared variance components cannot be silently omitted from Vp")
  disagree <- fit
  disagree$param_meta$component[1] <- "pen"
  check(all(is.na(subset(ge_parameters(ctx, disagree), parameter == "h2")$SE)),
        "Conflicting component metadata disables SE instead of mixing terms")
  native <- fit
  native$param_cov <- NULL
  native$native_se <- lapply(components, function(x) matrix(0.25, 4, 4))
  native_parameters <- ge_parameters(ctx, native)
  check(subset(native_parameters, parameter == "pen_variance" & trait == "TNB")$SE == 0.25 &&
          all(is.na(subset(native_parameters, parameter == "h2")$SE)),
        "Named native component SE remains available without fabricating derived covariance")

  zero_ctx <- list(nt = 1L, traits = ctx$traits[1, , drop = FALSE],
                   random = list(term("additive", "additive", 1L)),
                   constraints = list(additive = matrix(FALSE, 1, 1), residual = matrix(FALSE, 1, 1)))
  zero_fit <- list(components = list(additive = matrix(0, 1, 1), residual = matrix(0, 1, 1)),
                   param_meta = data.frame(role = c("additive", "residual"), trait1 = 1L, trait2 = 1L),
                   param_cov = diag(2), converged = TRUE)
  check(is.na(subset(ge_parameters(zero_ctx, zero_fit), parameter == "h2")$estimate),
        "Zero total variance produces undefined heritability rather than division artifacts")

  tiny_bad <- fit
  tiny_bad$param_cov <- diag(nrow(meta)) * 1e-12
  tiny_bad$param_cov[1, 2] <- tiny_bad$param_cov[2, 1] <- 2e-12
  check(all(is.na(subset(ge_parameters(ctx, tiny_bad), parameter == "h2")$SE)),
        "Small absolute scale cannot hide an indefinite sampling covariance")
  tiny_component <- fit
  tiny_component$components$birth_litter <- mat(1:2, c(1e-12, 1e-12), 2e-12)
  check(inherits(try(ge_parameters(ctx, tiny_component), silent = TRUE), "try-error"),
        "Small absolute scale cannot hide an indefinite variance component")
  units <- c(1e10, 1e-8, 1, 1)
  check(ge_result_psd(components$additive * outer(units, units)),
        "Valid covariance remains valid with strongly different trait units")

  named_iid <- ctx
  named_iid$random <- list(term("additive", "additive", 1:4), term("pe", "iid", 1:4))
  named_iid$constraints <- list(additive = ctx$constraints$additive, pe = matrix(FALSE, 4, 4), residual = ctx$constraints$residual)
  named_fit <- fit
  named_fit$components <- list(additive = components$additive, pe = mat(1:4, rep(1, 4)), residual = components$residual)
  named_fit$param_cov <- NULL
  check(!any(ge_parameters(named_iid, named_fit)$parameter == "repeatability"), "Name pe does not turn an iid term into PE")

  ctx$data <- data.frame(id = c(1L, 2L, 3L, 4L, 1L, 1L, 2L, 2L, 10L),
                         event_year = c(2020L, 2020L, 2020L, 2020L, 2021L, 2022L, 2021L, 2023L, 2020L),
                         birth_year = c(rep(2018L, 8), 2019L),
                         y1 = c(10, 11, 12, 13, rep(NA_real_, 5)),
                         y2 = c(1, 2, 3, NA, rep(NA_real_, 5)),
                         y3 = c(rep(NA_real_, 4), 13, 14, 15, 16, 20),
                         y4 = c(rep(NA_real_, 4), 5, 6, 7, 8, 9))
  ctx$trait_rows <- matrix(FALSE, nrow(ctx$data), 4)
  ctx$trait_rows[1:4, 1:2] <- TRUE
  ctx$trait_rows[5:9, 3:4] <- TRUE
  ctx$raw_data <- ctx$data
  ebv <- expand.grid(id = c(1:4, 10), trait = 1:4)
  ebv$EBV <- ebv$id + ebv$trait / 10
  ebv$original_id <- sprintf("%05d", ebv$id)
  ebv$SEP <- 0.2
  ebv$F <- 0
  ebv$has_record <- TRUE
  annual <- ge_annual(ctx, ebv)
  growth <- subset(annual$annual, trait == "age")
  check(nrow(growth) == 1L && growth$year == 2020 && growth$n_ebv == 4L,
        "Growth EBV grouped only in true growth events")
  backfat <- subset(annual$annual, trait == "bf")
  check(backfat$n_ebv == 4L && backfat$n_phenotyped == 3L,
        "Unmeasured backfat on a true growth event retains EBV cohort membership")
  check(!10L %in% ctx$data$id[ctx$trait_rows[, 1]] && abs(growth$mean_EBV - 2.6) < 1e-12,
        "Reproduction-only animal cannot contaminate growth annual EBV")
  check(subset(annual$annual, trait == "TNB" & year == 2021)$n_ebv == 2L,
        "Repeated-trait cohort follows its own event rows")
  bad_membership <- ctx
  bad_membership$trait_rows[1, 1] <- FALSE
  check(inherits(try(ge_annual(bad_membership, ebv), silent = TRUE), "try-error"),
        "Observed value outside event membership rejected")

  export <- ge_export_results(ctx, fit, ebv, file.path(output, "example_export"))
  check(file.exists(export$workbook) && file.exists(export$ebv_csv), "Generic workbook and EBV export complete")

  five <- list(nt = 5L, traits = data.frame(name = letters[1:5], label = LETTERS[1:5], unit = "unit"))
  five_annual <- expand.grid(year = 2020:2022, trait = letters[1:5], stringsAsFactors = FALSE)
  five_annual$year_basis <- "event_year"
  five_annual$mean_phenotype <- seq_len(nrow(five_annual))
  five_annual$mean_EBV <- five_annual$mean_phenotype / 10
  pages <- ge_result_trend_figures(five, five_annual, output)
  check(length(pages) == 4L && all(file.exists(pages)), "Five-trait trend plots paginate without oversized canvases")

  invisible(NULL)
}
