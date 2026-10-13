ge_result_or <- function(value, fallback) {
  if (is.null(value)) fallback else value
}

ge_result_components <- function(ctx, fit) {
  supplied <- fit$components

  if (is.null(supplied$additive) || is.null(supplied$residual)) {
    stop("Final additive and residual matrices are required.")
  }

  answer <- list(additive = as.matrix(supplied$additive))

  if (!is.null(ctx$shared_role)) {
    shared <- supplied$shared

    if (is.null(shared)) {
      shared <- supplied[[ctx$shared_role]]
    }

    if (is.null(shared)) {
      stop("The final shared-effect matrix is missing.")
    }

    answer[[ctx$shared_role]] <- as.matrix(shared)
  }

  answer$residual <- as.matrix(supplied$residual)

  for (role in names(answer)) {
    value <- answer[[role]]

    if (!identical(dim(value), c(as.integer(ctx$nt), as.integer(ctx$nt)))) {
      stop("A final component has the wrong matrix dimensions: ", role)
    }

    if (!is.numeric(value) || any(!is.finite(value)) ||
        max(abs(value - t(value))) > 1e-8 * max(1, abs(value))) {
      stop("A final component is non-finite or asymmetric: ", role)
    }

    if (min(eigen(value, symmetric = TRUE, only.values = TRUE)$values) <
        -1e-8 * max(1, abs(value))) {
      stop("A final component is not positive semidefinite: ", role)
    }
  }

  answer
}

ge_result_covariance <- function(ctx, fit, components) {
  meta <- fit$param_meta
  covariance <- fit$param_cov
  missing <- list(meta = NULL, covariance = NULL,
                  note = "Full sampling covariance is unavailable; derived SE was not calculated.")

  if (is.null(meta) || is.null(covariance)) {
    return(missing)
  }

  required <- c("role", "trait1", "trait2")

  if (!is.data.frame(meta) || !all(required %in% names(meta))) {
    missing$note <- "Sampling parameter metadata is incomplete; derived SE was not calculated."
    return(missing)
  }

  meta <- meta[required]
  meta$role <- as.character(meta$role)
  meta$role[meta$role == "shared"] <- ge_result_or(ctx$shared_role, "shared")
  meta$trait1 <- as.integer(meta$trait1)
  meta$trait2 <- as.integer(meta$trait2)
  covariance <- as.matrix(covariance)

  valid_meta <- all(meta$role %in% names(components)) &&
    all(!is.na(meta$trait1) & !is.na(meta$trait2)) &&
    all(meta$trait1 >= 1 & meta$trait1 <= ctx$nt) &&
    all(meta$trait2 >= 1 & meta$trait2 <= ctx$nt)

  keys <- paste(meta$role, pmin(meta$trait1, meta$trait2),
                pmax(meta$trait1, meta$trait2), sep = ":")

  if (!valid_meta || !nrow(meta) || anyDuplicated(keys) ||
      !identical(dim(covariance), c(nrow(meta), nrow(meta))) || !is.numeric(covariance) ||
      any(!is.finite(covariance))) {
    missing$note <- "Sampling covariance and parameter order could not be verified; derived SE was not calculated."
    return(missing)
  }

  for (role in names(components)) {
    constraint_role <- if (role %in% c("litter", "pe")) "shared" else role
    constraints <- ctx$constraints[[constraint_role]]

    for (first in seq_len(ctx$nt)) {
      for (second in first:ctx$nt) {
        fixed <- !is.null(constraints) && isTRUE(constraints[first, second])
        key <- paste(role, first, second, sep = ":")

        if ((!fixed && !key %in% keys) || (fixed && key %in% keys)) {
          missing$note <- "Sampling covariance does not match the estimated-parameter set; derived SE was not calculated."
          return(missing)
        }
      }
    }
  }

  scale <- max(1, abs(covariance))

  if (max(abs(covariance - t(covariance))) > 1e-8 * scale ||
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * scale) {
    missing$note <- "Sampling covariance is not a valid symmetric positive-semidefinite matrix; derived SE was not calculated."
    return(missing)
  }

  list(meta = meta, covariance = covariance, note = "")
}

ge_result_delta <- function(gradient, sampling) {
  if (is.null(sampling$covariance)) {
    return(NA_real_)
  }

  variance <- drop(crossprod(gradient, sampling$covariance %*% gradient))
  tolerance <- 1e-10 * max(1, abs(sampling$covariance)) * sum(abs(gradient))^2

  if (!is.finite(variance) || variance < -tolerance) {
    return(NA_real_)
  }

  sqrt(max(0, variance))
}

ge_result_gradient <- function(kind, components, meta, trait1, trait2 = trait1, role = NULL) {
  if (is.null(meta)) {
    return(numeric())
  }

  diagonal <- meta$trait1 == trait1 & meta$trait2 == trait1
  genetic <- diagonal & meta$role == "additive"
  total <- sum(vapply(components, function(x) x[trait1, trait1], numeric(1)))
  additive <- components$additive[trait1, trait1]
  gradient <- rep(0, nrow(meta))

  if (kind == "component") {
    same_pair <- pmin(meta$trait1, meta$trait2) == min(trait1, trait2) &
      pmax(meta$trait1, meta$trait2) == max(trait1, trait2)

    gradient[meta$role == role & same_pair] <- 1
  }

  if (kind == "Vp") {
    gradient[diagonal] <- 1
  }

  if (kind == "h2") {
    gradient[diagonal] <- -additive / total^2
    gradient[genetic] <- gradient[genetic] + 1 / total
  }

  if (kind == "repeatability") {
    numerator <- additive + components$pe[trait1, trait1]
    gradient[diagonal] <- -numerator / total^2
    numerator_parameters <- diagonal & meta$role %in% c("additive", "pe")
    gradient[numerator_parameters] <- gradient[numerator_parameters] + 1 / total
  }

  if (kind == "rg") {
    first <- components$additive[trait1, trait1]
    second <- components$additive[trait2, trait2]
    cross <- components$additive[trait1, trait2]
    correlation <- cross / sqrt(first * second)

    gradient[meta$role == "additive" & meta$trait1 == trait1 & meta$trait2 == trait1] <-
      -correlation / (2 * first)

    gradient[meta$role == "additive" & meta$trait1 == trait2 & meta$trait2 == trait2] <-
      -correlation / (2 * second)

    pair <- pmin(meta$trait1, meta$trait2) == min(trait1, trait2) &
      pmax(meta$trait1, meta$trait2) == max(trait1, trait2)

    gradient[meta$role == "additive" & pair] <- 1 / sqrt(first * second)
  }

  gradient
}

ge_result_native_se <- function(fit, role, trait1, trait2) {
  value <- fit$native_se

  if (!is.list(value)) {
    return(NA_real_)
  }

  component <- value[[role]]

  if (is.null(component) && role %in% c("litter", "pe")) {
    component <- value$shared
  }

  if (is.null(component)) {
    return(NA_real_)
  }

  component <- as.matrix(component)

  if (nrow(component) < trait1 || ncol(component) < trait2) {
    return(NA_real_)
  }

  answer <- suppressWarnings(as.numeric(component[trait1, trait2]))
  if (is.finite(answer) && answer >= 0) answer else NA_real_
}

ge_parameters <- function(ctx, fit) {
  components <- ge_result_components(ctx, fit)
  sampling <- ge_result_covariance(ctx, fit, components)
  rows <- list()

  add_row <- function(parameter, first, second = first, estimate,
                      kind, role = NULL, note = "") {
    gradient <- ge_result_gradient(kind, components, sampling$meta, first, second, role)
    standard_error <- if (is.finite(estimate)) ge_result_delta(gradient, sampling) else NA_real_
    method <- if (is.finite(standard_error)) "delta" else "unavailable"

    if (kind == "component" && !is.finite(standard_error)) {
      standard_error <- ge_result_native_se(fit, role, first, second)

      if (is.finite(standard_error)) {
        method <- "software"
      }
    }

    if (kind == "component" && method == "delta") {
      method <- "sampling covariance"
    }

    if (!is.finite(standard_error) && is.finite(estimate)) {
      note <- paste(note, sampling$note)
    }

    dimensionless <- kind %in% c("h2", "repeatability", "rg")
    unit <- if (dimensionless) "ratio" else if (first == second) {
      paste0(ctx$traits$unit[first], "^2")
    } else {
      paste(ctx$traits$unit[first], ctx$traits$unit[second], sep = " * ")
    }

    trait <- if (first == second) as.character(ctx$traits$name[first]) else {
      paste(ctx$traits$name[c(first, second)], collapse = " : ")
    }

    if (!isTRUE(fit$converged)) {
      note <- paste(note, "Estimation did not converge; values are provisional.")
    }

    rows[[length(rows) + 1L]] <<- data.frame(
      parameter = parameter, trait = trait, estimate = estimate,
      SE = standard_error, unit = unit, SE_method = method,
      note = trimws(note), stringsAsFactors = FALSE
    )
  }

  for (trait in seq_len(ctx$nt)) {
    for (role in names(components)) {
      estimate <- components[[role]][trait, trait]
      note <- if (estimate <= 1e-10) "Boundary estimate; asymptotic SE requires caution." else ""

      add_row(paste0(role, "_variance"), trait, estimate = estimate,
              kind = "component", role = role, note = note)
    }

    total <- sum(vapply(components, function(x) x[trait, trait], numeric(1)))
    boundary <- any(vapply(components, function(x) x[trait, trait] <= 1e-10, logical(1)))
    note <- if (boundary) "A variance component is on the boundary; delta SE is an asymptotic approximation." else ""

    add_row("Vp", trait, estimate = total, kind = "Vp", note = note)

    if (total > 0) {
      add_row("h2", trait, estimate = components$additive[trait, trait] / total,
              kind = "h2", note = note)

      if ("pe" %in% names(components)) {
        add_row("repeatability", trait,
                estimate = (components$additive[trait, trait] + components$pe[trait, trait]) / total,
                kind = "repeatability", note = note)
      }
    } else {
      add_row("h2", trait, estimate = NA_real_, kind = "h2",
              note = "Total variance is zero; heritability is undefined.")

      if ("pe" %in% names(components)) {
        add_row("repeatability", trait, estimate = NA_real_, kind = "repeatability",
                note = "Total variance is zero; repeatability is undefined.")
      }
    }
  }

  if (ctx$nt > 1) {
    pairs <- combn(seq_len(ctx$nt), 2)

    for (pair in seq_len(ncol(pairs))) {
      first <- pairs[1, pair]
      second <- pairs[2, pair]

      for (role in names(components)) {
        constraint_role <- if (role %in% c("litter", "pe")) "shared" else role
        constraint <- ctx$constraints[[constraint_role]]
        fixed <- !is.null(constraint) && isTRUE(constraint[first, second])
        note <- if (fixed) "Fixed zero covariance; not an estimated parameter." else ""

        add_row(paste0(role, "_covariance"), first, second,
                estimate = components[[role]][first, second], kind = "component",
                role = role, note = note)

        if (fixed) {
          rows[[length(rows)]]$SE <- NA_real_
          rows[[length(rows)]]$SE_method <- "fixed"
        }
      }

      first_variance <- components$additive[first, first]
      second_variance <- components$additive[second, second]
      valid <- first_variance > 0 && second_variance > 0
      correlation <- if (valid) {
        components$additive[first, second] / sqrt(first_variance * second_variance)
      } else NA_real_

      note <- if (!valid) "Additive variance is zero; genetic correlation is undefined." else {
        if (abs(correlation) >= 1 - 1e-8) "Boundary correlation; delta SE requires caution." else ""
      }

      add_row("rg", first, second, estimate = correlation, kind = "rg", note = note)
    }
  }

  do.call(rbind, rows)
}

ge_result_ebv <- function(ctx, ebv) {
  answer <- as.data.frame(ebv, stringsAsFactors = FALSE)
  required <- c("id", "trait", "EBV")

  if (!all(required %in% names(answer))) {
    stop("EBV output must contain id, trait and EBV.")
  }

  if (anyNA(answer$id) || anyNA(answer$trait) ||
      any(!answer$trait %in% seq_len(ctx$nt)) ||
      anyDuplicated(answer[c("id", "trait")])) {
    stop("EBV output has missing, duplicate or invalid animal/trait identifiers.")
  }

  answer
}

ge_annual <- function(ctx, ebv) {
  ebv <- ge_result_ebv(ctx, ebv)
  bases <- ge_result_or(ctx$trend_basis, "event_year")
  bases <- unique(as.character(bases))

  if (!all(bases %in% c("event_year", "birth_year"))) {
    stop("trend_basis must be event_year and/or birth_year.")
  }

  rows <- list()
  warnings <- character()

  for (basis in bases) {
    if (!basis %in% names(ctx$data)) {
      warnings <- c(warnings, paste("Annual summary unavailable:", basis, "is missing."))
      next
    }

    raw_year <- as.character(ctx$data[[basis]])
    valid_year <- !is.na(raw_year) & grepl("^[0-9]{4}$", raw_year)
    years <- suppressWarnings(as.integer(raw_year))
    years[!valid_year] <- NA_integer_

    if (any(!valid_year)) {
      warnings <- c(warnings, paste(sum(!valid_year), "records omitted from", basis, "summary because the year is unavailable."))
    }

    if (basis == "birth_year") {
      id_year <- unique(data.frame(id = ctx$data$id[valid_year], year = years[valid_year]))

      if (anyDuplicated(id_year$id)) {
        stop("An animal has conflicting birth years.")
      }
    }

    for (year in sort(unique(years[valid_year]))) {
      subset <- ctx$data[!is.na(years) & years == year, , drop = FALSE]
      ids <- unique(subset$id)

      for (trait in seq_len(ctx$nt)) {
        phenotype <- subset[[paste0("y", trait)]]
        observed <- is.finite(phenotype)
        genetic <- ebv$EBV[ebv$trait == trait & ebv$id %in% ids]
        genetic <- genetic[is.finite(genetic)]

        rows[[length(rows) + 1L]] <- data.frame(
          year_basis = basis, year = year, trait = as.character(ctx$traits$name[trait]),
          n_records = sum(observed), n_animals = length(ids),
          n_phenotyped = length(unique(subset$id[observed])), n_ebv = length(genetic),
          mean_phenotype = if (any(observed)) mean(phenotype[observed]) else NA_real_,
          sd_phenotype = if (sum(observed) > 1) sd(phenotype[observed]) else NA_real_,
          mean_EBV = if (length(genetic)) mean(genetic) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  annual <- if (length(rows)) do.call(rbind, rows) else data.frame(
    year_basis = character(), year = integer(), trait = character(),
    n_records = integer(), n_animals = integer(), n_phenotyped = integer(), n_ebv = integer(),
    mean_phenotype = numeric(), sd_phenotype = numeric(), mean_EBV = numeric()
  )

  trends <- list()

  for (basis in unique(annual$year_basis)) {
    for (trait in seq_len(ctx$nt)) {
      subset <- annual[annual$year_basis == basis & annual$trait == ctx$traits$name[trait], , drop = FALSE]

      for (type in c("phenotypic", "genetic")) {
        value <- subset[[if (type == "phenotypic") "mean_phenotype" else "mean_EBV"]]
        keep <- is.finite(value) & is.finite(subset$year)
        n_years <- sum(keep)
        slope <- standard_error <- r_squared <- NA_real_

        if (n_years >= 3) {
          model <- lm(value[keep] ~ subset$year[keep])
          model_summary <- summary(model)
          slope <- unname(coef(model)[2])
          standard_error <- unname(model_summary$coefficients[2, 2])
          r_squared <- model_summary$r.squared

          if (!is.finite(r_squared)) {
            r_squared <- NA_real_
          }
        }

        trends[[length(trends) + 1L]] <- data.frame(
          trait = as.character(ctx$traits$name[trait]), trend = type,
          year_basis = basis, n_years = n_years, change_per_year = slope,
          SE = standard_error, R_squared = r_squared,
          unit = paste0(ctx$traits$unit[trait], "/year"),
          note = if (n_years < 3) "Fewer than three effective years; no linear trend fitted." else
            "Unweighted regression of annual means; EBV counted once per animal per year.",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  trend_table <- if (length(trends)) do.call(rbind, trends) else data.frame(
    trait = character(), trend = character(), year_basis = character(), n_years = integer(),
    change_per_year = numeric(), SE = numeric(), R_squared = numeric(), unit = character(), note = character()
  )

  list(annual = annual, trends = trend_table, warnings = warnings)
}

ge_result_reliability <- function(ctx, fit, ebv) {
  answer <- ge_result_ebv(ctx, ebv)
  answer$reliability <- NA_real_
  answer$reliability_status <- "Verified SEP and animal relationship diagonal are required."

  if (!all(c("SEP", "F") %in% names(answer))) {
    return(answer)
  }

  additive <- ge_result_components(ctx, fit)$additive

  for (row in seq_len(nrow(answer))) {
    trait <- answer$trait[row]
    sep <- answer$SEP[row]
    inbreeding <- answer$F[row]
    variance <- additive[trait, trait]

    if (!is.finite(sep) || sep < 0 || !is.finite(inbreeding) ||
        inbreeding < 0 || inbreeding > 1 || variance <= 0) {
      next
    }

    individual_variance <- (1 + inbreeding) * variance
    reliability <- 1 - sep^2 / individual_variance

    if (reliability < 0 || reliability > 1) {
      answer$reliability_status[row] <- "PEV is inconsistent with the verified model variance; value was not clipped."
      next
    }

    answer$reliability[row] <- reliability
    answer$reliability_status[row] <- "Verified: 1 - SEP^2 / ((1 + F) * additive variance)."
  }

  answer$trait_name <- as.character(ctx$traits$name[answer$trait])
  answer
}

ge_result_write_csv <- function(value, file) {
  connection <- file(file, open = "wb")
  on.exit(close(connection))
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), connection)
  lines <- capture.output(write.table(value, sep = ",", row.names = FALSE,
                                     col.names = TRUE, na = "", qmethod = "double"))
  writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\r\n"), "\r\n")), connection)
}

ge_result_sheet <- function(workbook, name, value, digits = 3) {
  if (is.null(value)) {
    return(invisible(NULL))
  }

  if (!is.data.frame(value)) {
    value <- as.data.frame(value, stringsAsFactors = FALSE)
  }

  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, value, keepNA = FALSE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)

  header <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#D35C2D",
                                 textDecoration = "bold", halign = "center", wrapText = TRUE)

  if (ncol(value)) {
    openxlsx::addStyle(workbook, name, header, rows = 1, cols = seq_len(ncol(value)), gridExpand = TRUE)
    openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(value)), widths = "auto")
  }

  if (nrow(value)) {
    for (column in seq_len(ncol(value))) {
      if (!is.numeric(value[[column]])) {
        next
      }

      label <- names(value)[column]
      count <- is.integer(value[[column]]) ||
        grepl("(^n(_|$)|count|records|animals|sows|year$|rounds|^trait$|^id$|number)", label,
              ignore.case = TRUE) ||
        grepl("\u8bb0\u5f55\u6570|\u4e2a\u4f53\u6570|\u6bcd\u732a\u6570|\u5e74\u4efd", label)
      format <- if (count) "0" else paste0("0.", paste(rep("0", digits), collapse = ""))
      style <- openxlsx::createStyle(numFmt = format)
      openxlsx::addStyle(workbook, name, style, rows = seq_len(nrow(value)) + 1L,
                        cols = column, gridExpand = TRUE, stack = TRUE)
    }
  }

  invisible(NULL)
}

ge_result_descriptive <- function(ctx) {
  answer <- ctx$descriptive

  if (!is.null(answer)) {
    return(as.data.frame(answer, stringsAsFactors = FALSE))
  }

  rows <- lapply(seq_len(ctx$nt), function(trait) {
    value <- ctx$data[[paste0("y", trait)]]
    value <- value[is.finite(value)]
    average <- if (length(value)) mean(value) else NA_real_
    deviation <- if (length(value) > 1) sd(value) else NA_real_

    data.frame(trait = ctx$traits$name[trait], stage = "post_QC", n = length(value),
               mean = average, sd = deviation,
               median = if (length(value)) median(value) else NA_real_,
               min = if (length(value)) min(value) else NA_real_,
               max = if (length(value)) max(value) else NA_real_,
               CV_percent = if (is.finite(average) && average > 0) 100 * deviation / average else NA_real_)
  })

  do.call(rbind, rows)
}

ge_result_png <- function(file, width, height, draw) {
  grDevices::png(file, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.5, 4.5, 3, 1), family = "sans")
  draw()
  invisible(file)
}

ge_result_histograms <- function(ctx, directory) {
  before <- ge_result_or(ctx$raw_data, ctx$qc_before)
  paths <- character()

  for (trait in seq_len(ctx$nt)) {
    column <- paste0("y", trait)
    after <- ctx$data[[column]]
    after <- after[is.finite(after)]
    original <- if (!is.null(before) && column %in% names(before)) before[[column]] else numeric()
    original <- original[is.finite(original)]
    file <- file.path(directory, paste0("histogram_", trait, ".png"))
    values <- c(original, after)

    ge_result_png(file, if (length(original)) 1600 else 900, 700, function() {
      if (length(original)) {
        graphics::par(mfrow = c(1, 2))
      }

      groups <- if (length(original)) list(original, after) else list(after)
      labels <- if (length(original)) c("Before QC", "After QC") else "After QC"
      breaks <- if (length(values) && diff(range(values)) > 0) {
        pretty(range(values), n = 20)
      } else if (length(values)) values[1] + c(-0.5, 0.5) else c(0, 1)
      limits <- if (length(values)) range(breaks) else c(0, 1)
      frequencies <- lapply(groups, function(x) if (length(x)) hist(x, breaks = breaks, plot = FALSE)$counts else 0)
      ymax <- max(1, unlist(frequencies)) * 1.1

      for (index in seq_along(groups)) {
        title <- paste0(ctx$traits$label[trait], ": ", labels[index], "\nn = ", length(groups[[index]]))

        if (length(groups[[index]])) {
          graphics::hist(groups[[index]], breaks = breaks, xlim = limits, ylim = c(0, ymax),
                         col = if (labels[index] == "Before QC") "#C8C8C8" else "#4A83AE",
                         border = "white", main = title,
                         xlab = paste0(ctx$traits$label[trait], " (", ctx$traits$unit[trait], ")"),
                         ylab = "Records")
        } else {
          graphics::plot.new()
          graphics::title(main = title)
          graphics::text(0.5, 0.5, "No valid observations")
        }
      }
    })

    paths <- c(paths, file)
  }

  paths
}

ge_result_trend_figures <- function(ctx, annual, directory) {
  paths <- character()

  for (basis in unique(annual$year_basis)) {
    for (type in c("phenotypic", "genetic")) {
      file <- file.path(directory, paste0(type, "_", basis, ".png"))
      column <- if (type == "phenotypic") "mean_phenotype" else "mean_EBV"
      title <- if (type == "phenotypic") "\u8868\u578b\u8fdb\u5c55" else "\u9057\u4f20\u8fdb\u5c55"

      ge_result_png(file, 1400, 700 * ctx$nt, function() {
        graphics::par(mfrow = c(ctx$nt, 1))

        for (trait in seq_len(ctx$nt)) {
          subset <- annual[annual$year_basis == basis & annual$trait == ctx$traits$name[trait], , drop = FALSE]
          keep <- is.finite(subset[[column]]) & is.finite(subset$year)
          value <- subset[[column]][keep]
          year <- subset$year[keep]

          if (!length(value)) {
            graphics::plot.new()
            graphics::title(main = title)
            graphics::text(0.5, 0.5, "No valid annual means")
            next
          }

          xlab <- if (basis == "birth_year") "Birth year" else "Production year"
          graphics::plot(year, value, type = "o", pch = 19, lwd = 2,
                         col = "#D35C2D", main = title, xlab = xlab,
                         ylab = paste0(ctx$traits$label[trait], " (", ctx$traits$unit[trait], ")"), xaxt = "n")
          graphics::axis(1, at = year, labels = year)

          if (length(value) >= 3) {
            graphics::abline(lm(value ~ year), col = "#595959", lty = 2, lwd = 1.5)
          }
        }
      })

      paths <- c(paths, file)
    }
  }

  paths
}

ge_export_results <- function(ctx, fit, ebv, directory) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx is required to export the result workbook.")
  }

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  parameters <- ge_parameters(ctx, fit)
  annual <- ge_annual(ctx, ebv)
  ebv <- ge_result_reliability(ctx, fit, ebv)
  descriptive <- ge_result_descriptive(ctx)

  settings <- data.frame(
    item = c("analysis", "species", "model_family", "software", "version", "trait_order",
             "shared_effect", "trend_basis", "converged", "iterations", "minus2logL", "model"),
    value = c(as.character(ctx$name), as.character(ctx$species), as.character(ctx$family),
              as.character(fit$engine), paste(fit$version, collapse = "; "),
              paste(ctx$traits$name, collapse = ", "), ge_result_or(ctx$shared_role, "none"),
              paste(ge_result_or(ctx$trend_basis, "event_year"), collapse = ", "),
              as.character(fit$converged), as.character(ge_result_or(fit$rounds, NA_integer_)),
              as.character(ge_result_or(fit$minus2logL, NA_real_)),
              paste(ge_result_or(ctx$model_description, "See retained model configuration and parameter cards."), collapse = "\n")),
    stringsAsFactors = FALSE
  )

  warnings <- c(as.character(fit$warnings), annual$warnings,
                if (is.null(ctx$raw_data) && is.null(ctx$qc_before))
                  "Pre-QC histogram data unavailable; only retained distributions were plotted.")

  diagnostics <- data.frame(message = unique(warnings[!is.na(warnings) & nzchar(warnings)]), stringsAsFactors = FALSE)
  workbook <- openxlsx::createWorkbook()
  ge_result_sheet(workbook, "Analysis", settings)
  ge_result_sheet(workbook, "Descriptive", descriptive, 2)
  ge_result_sheet(workbook, "QC", ctx$qc_audit, 2)
  ge_result_sheet(workbook, "Input audit", ctx$input_audit, 2)
  ge_result_sheet(workbook, "Genetic parameters", parameters, 3)
  ge_result_sheet(workbook, "Annual means", annual$annual, 2)
  ge_result_sheet(workbook, "Annual trends", annual$trends, 3)
  ge_result_sheet(workbook, "Diagnostics", diagnostics)

  if (nrow(ebv) <= 1048575L) {
    ge_result_sheet(workbook, "EBV", ebv, 3)
  }

  xlsx <- file.path(directory, "results.xlsx")
  csv <- file.path(directory, "EBV.csv")
  openxlsx::saveWorkbook(workbook, xlsx, overwrite = TRUE)
  ge_result_write_csv(ebv, csv)
  histograms <- ge_result_histograms(ctx, directory)
  trend_figures <- ge_result_trend_figures(ctx, annual$annual, directory)

  list(workbook = xlsx, ebv_csv = csv, histograms = histograms,
       trend_figures = trend_figures, warnings = warnings,
       parameters = parameters, annual = annual$annual, trends = annual$trends)
}
