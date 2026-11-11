script_path <- if (sys.nframe() > 0L) {
  sys.frame(1)$ofile
} else {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  sub("^--file=", "", file_arg[1L])
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "TDA.R"))

main_scenarios <- function() {
  position <- expand.grid(
    design = c("EQ", "NI"),
    mu_ratio = c(0, 0.5, 0.9, 1.1),
    tau_ratio = c(0.05, 0.10, 0.20),
    ratio = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  position$module <- "position"
  position$margin_pattern <- "constant"

  margins <- expand.grid(
    design = c("EQ", "NI"),
    mu_ratio = c(0.9, 1.1),
    tau_ratio = 0.10,
    margin_pattern = c("shuffled", "increasing", "decreasing"),
    ratio = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  margins$module <- "margins"

  scenarios <- rbind(position, margins[, names(position)])
  scenarios$stage <- "primary"
  scenarios$precision <- "main"
  scenarios$distribution <- "normal"
  scenarios$k <- 50L
  scenarios$id <- seq_len(nrow(scenarios))
  scenarios$delta_ref <- 1
  scenarios$mu <- scenarios$mu_ratio * scenarios$delta_ref
  scenarios$tau <- scenarios$tau_ratio * scenarios$delta_ref
  scenarios$trial_alpha <- 0.05
  scenarios$test_alpha <- 0.10
  scenarios$reps <- ifelse(scenarios$ratio == 1, 1000L, 500L)
  scenarios$seed_base <- 721000000L + 100000L * scenarios$id
  scenarios
}

distribution_scenarios <- function() {
  main <- expand.grid(
    design = c("EQ", "NI"), k = c(20L, 50L),
    truth = c("normal", "t4", "skew"), ratio = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  main$family <- "main"
  main$mu <- 0
  main$se_pattern <- "variable"
  main$slope <- 0

  boundary <- expand.grid(
    design = c("EQ", "NI"), k = c(20L, 50L),
    truth = c("point", "near_merge", "variance_boundary", "interior"),
    stringsAsFactors = FALSE
  )
  boundary$ratio <- 1
  boundary$family <- "boundary"
  boundary$mu <- 0
  boundary$se_pattern <- "variable"
  boundary$slope <- 0

  equal_se <- expand.grid(
    design = c("EQ", "NI"), ratio = c(1, 5),
    stringsAsFactors = FALSE
  )
  equal_se$k <- 50L
  equal_se$truth <- "normal"
  equal_se$family <- "constant_se"
  equal_se$mu <- 0
  equal_se$se_pattern <- "constant"
  equal_se$slope <- 0

  margin <- data.frame(
    design = "NI", k = 50L, truth = "normal", ratio = c(1, 2),
    family = "ni_margin", mu = 1, se_pattern = "variable", slope = 0
  )
  dependent <- data.frame(
    design = c("EQ", "NI"), k = 50L, truth = "normal", ratio = 1,
    family = "design_dependence", mu = 0, se_pattern = "variable", slope = 1
  )

  columns <- c("family", "design", "k", "truth", "ratio", "mu",
               "se_pattern", "slope")
  scenarios <- do.call(rbind, lapply(
    list(main, boundary, equal_se, margin, dependent),
    function(x) x[, columns]
  ))
  scenarios$id <- seq_len(nrow(scenarios))
  scenarios$reps <- ifelse(
    scenarios$ratio == 1,
    ifelse(scenarios$family == "main", 1000L, 500L),
    300L
  )
  scenarios
}

tda_scenarios <- function() {
  k50 <- main_scenarios()
  k30 <- k50
  k30$k <- 30L
  k30$seed_base <- 941000000L + 100000L * k30$id
  distribution <- distribution_scenarios()

  k50$family <- k30$family <- "primary"
  k50$truth <- k30$truth <- "normal"
  k50$archive <- "primary_k50"
  k30$archive <- "primary_k30"
  k50$se_pattern <- k30$se_pattern <- ""
  k50$slope <- k30$slope <- NA_real_
  distribution$archive <- "distribution"
  distribution$seed_base <- 311000000L + 100000L * distribution$id
  distribution$trial_alpha <- 0.05

  columns <- c(
    "archive", "family", "design", "k", "id", "module", "mu_ratio",
    "tau_ratio", "delta_ref", "mu", "tau", "margin_pattern", "precision",
    "distribution", "truth", "ratio", "se_pattern", "slope", "trial_alpha",
    "seed_base", "reps"
  )
  fill_columns <- function(data) {
    for (name in setdiff(columns, names(data))) data[[name]] <- NA
    data[, columns]
  }
  scenarios <- rbind(
    fill_columns(k50), fill_columns(k30), fill_columns(distribution)
  )
  rownames(scenarios) <- NULL
  scenarios$cell_index <- seq_len(nrow(scenarios))
  scenarios$original_reps <- scenarios$reps
  scenarios$reps <- 1000L
  scenarios$test_alpha <- 0.10
  scenarios$work_nodes <- 24L
  scenarios
}

margin_permutations <- new.env(parent = emptyenv())

shuffled_margins <- function(values, s) {
  key <- as.character(length(s))
  if (!exists(key, envir = margin_permutations, inherits = FALSE)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)

    set.seed(721903L)
    permutations <- replicate(1000L, sample(values))
    correlations <- apply(permutations, 2L, function(x) cor(x, log(s)))
    selected <- permutations[, which.min(abs(correlations))]
    assign(key, selected, envir = margin_permutations)

    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
  get(key, envir = margin_permutations, inherits = FALSE)
}

main_design <- function(scenario) {
  ref <- scenario$delta_ref
  endpoints <- switch(
    scenario$precision,
    main = c(0.25, 0.45), high = c(0.10, 0.30), low = c(0.40, 0.65)
  )
  s_ratio <- exp(seq(log(endpoints[1]), log(endpoints[2]),
                     length.out = scenario$k))
  values <- seq(0.8, 1.2, length.out = scenario$k)
  margin_ratio <- switch(
    scenario$margin_pattern,
    constant = rep(1, scenario$k),
    increasing = values,
    decreasing = rev(values),
    shuffled = shuffled_margins(values, s_ratio)
  )

  delta <- ref * margin_ratio
  s <- ref * s_ratio
  critical <- qnorm(1 - scenario$trial_alpha)
  lower <- if (scenario$design == "EQ") {
    -delta + critical * s
  } else {
    rep(-Inf, scenario$k)
  }
  list(
    s = s, delta = delta, lower = lower, upper = delta - critical * s,
    mu = scenario$mu_ratio * ref, tau = scenario$tau_ratio * ref,
    delta_ref = ref
  )
}

distribution_design <- function(scenario) {
  s <- if (scenario$se_pattern == "constant") {
    rep(0.4, scenario$k)
  } else {
    exp(seq(log(0.20), log(0.65), length.out = scenario$k))
  }
  shift <- scenario$mu + scenario$slope * (log(s) - mean(log(s)))
  lower <- if (scenario$design == "EQ") {
    -1 + qnorm(0.95) * s
  } else {
    rep(-Inf, scenario$k)
  }
  list(s = s, lo = lower, hi = 1 - qnorm(0.95) * s, shift = shift)
}

tda_sim_design <- function(scenario) {
  if (scenario$family == "primary") {
    design <- main_design(scenario)
    design$lo <- design$lower
    design$hi <- design$upper
    design
  } else {
    distribution_design(scenario)
  }
}

generate_main <- function(scenario, design) {
  vapply(seq_len(scenario$k), function(i) {
    repeat {
      z <- switch(
        scenario$distribution,
        normal = rnorm(1),
        t4 = rt(1, 4) / sqrt(2),
        skew = (exp(rnorm(1)) - exp(0.5)) / sqrt(exp(1) * (exp(1) - 1))
      )
      theta <- design$mu + design$tau * z
      y <- theta + rnorm(1, 0, design$s[i])
      favorable <- design$lower[i] < y && y < design$upper[i]
      publication_probability <- if (favorable) 1 else 1 / scenario$ratio
      if (runif(1) < publication_probability) return(y)
    }
  }, 0.0)
}

generate_distribution <- function(scenario, design) {
  vapply(seq_along(design$s), function(i) {
    repeat {
      theta <- switch(
        scenario$truth,
        normal = rnorm(1, 0, 0.4),
        t4 = 0.4 * rt(1, 4) / sqrt(2),
        skew = 0.4 * (exp(rnorm(1)) - exp(0.5)) /
          sqrt(exp(1) * (exp(1) - 1)),
        point = 0,
        near_merge = rnorm(1, if (runif(1) < 0.5) -0.05 else 0.05,
                           sqrt(0.1575)),
        variance_boundary = {
          first_component <- runif(1) < 0.35
          rnorm(1, if (first_component) -0.4 else 0.5,
                sqrt(if (first_component) 0 else 0.30))
        },
        interior = {
          first_component <- runif(1) < 0.35
          rnorm(1, if (first_component) -0.4 else 0.5,
                sqrt(if (first_component) 0.08 else 0.30))
        }
      )
      y <- theta + design$shift[i] + rnorm(1, 0, design$s[i])
      favorable <- y > design$lo[i] && y < design$hi[i]
      publication_probability <- if (favorable) 1 else 1 / scenario$ratio
      if (runif(1) < publication_probability) return(y)
    }
  }, 0.0)
}

tda_sim_generate <- function(scenario, replicate,
                             design = tda_sim_design(scenario)) {
  if (length(replicate) != 1L || !is.finite(replicate) || replicate < 1 ||
      replicate != as.integer(replicate) || replicate > 1000) {
    stop("replicate must be an integer between 1 and 1000")
  }
  seed <- if (replicate <= scenario$original_reps) {
    scenario$seed_base + 10L * replicate
  } else {
    1870000000L + 10000L * scenario$cell_index + replicate
  }
  set.seed(as.integer(seed), kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")

  y <- if (scenario$family == "primary") {
    generate_main(scenario, design)
  } else {
    generate_distribution(scenario, design)
  }
  list(y = y, seed = seed)
}

simulation_record <- function(scenario, replicate, sample, design, fit) {
  diagnostics <- fit$diagnostics
  data.frame(
    cell_index = scenario$cell_index,
    replicate = replicate,
    seed = sample$seed,
    p_tda = fit$p.value,
    p_egger = egger_p(sample$y, design$s, scenario$design),
    S = fit$S,
    V = fit$variance,
    status = fit$status,
    work_mu = fit$working_parameters[1L],
    work_tau2 = fit$working_parameters[2L],
    nonempty_regions = fit$nonempty_regions,
    informative_targets = fit$informative_targets,
    positive_pairs = fit$positive_pairs,
    fit_ok = diagnostics$fit_ok,
    numeric_ok = diagnostics$numeric_ok,
    fit_gap = diagnostics$fit_gap,
    cap_active = diagnostics$cap_active,
    max_load = diagnostics$max_load,
    residual = diagnostics$kkt_residual,
    refined = diagnostics$refined,
    unresolved = diagnostics$unresolved,
    row.names = NULL
  )
}

summarize_cell <- function(scenario, results) {
  summaries <- lapply(c("tda", "egger"), function(method) {
    p <- results[[paste0("p_", method)]]
    n <- sum(is.finite(p))
    rejection <- if (n) mean(p[is.finite(p)] <= 0.10) else NA_real_
    mcse <- if (n) 100 * sqrt(rejection * (1 - rejection) / n) else NA_real_

    data.frame(
      scenario,
      method = method,
      completed_replicates = nrow(results),
      defined_p_values = n,
      rejection_percent = 100 * rejection,
      mcse_percent = mcse,
      numerical_failures = if (method == "tda") {
        sum(results$status == "numerical_failure")
      } else 0L,
      no_information = if (method == "tda") {
        sum(results$status == "no_usable_cross_precision_information")
      } else 0L
    )
  })
  do.call(rbind, summaries)
}

tda_simulation <- function(cells = seq_len(276L), replicates = 1000L,
                           output_dir = "simulation_results", resume = TRUE) {
  scenarios <- tda_scenarios()
  if (!length(cells) || anyNA(cells) || any(cells != as.integer(cells)) ||
      any(!cells %in% scenarios$cell_index) || anyDuplicated(cells)) {
    stop("cells must contain unique scenario indices between 1 and 276")
  }
  if (length(replicates) != 1L || !is.finite(replicates) || replicates < 1L ||
      replicates > 1000L || replicates != as.integer(replicates)) {
    stop("replicates must be an integer between 1 and 1000")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  block_dir <- file.path(output_dir, "blocks")
  dir.create(block_dir, showWarnings = FALSE)
  selected <- scenarios[cells, ]
  selected$reps <- as.integer(replicates)
  write.csv(selected, file.path(output_dir, "scenarios.csv"), row.names = FALSE)

  sources <- file.path(script_dir, c("TDA.R", "simulation.R"))
  hashes <- tools::md5sum(sources)
  names(hashes) <- basename(names(hashes))
  summaries <- list()

  for (i in cells) {
    scenario <- scenarios[i, ]
    design <- tda_sim_design(scenario)
    rows <- list()
    initial <- numeric()

    for (start in seq.int(1L, replicates, by = 25L)) {
      end <- min(start + 24L, replicates)
      filename <- sprintf("cell_%03d_reps_%04d_%04d.rds", i, start, end)
      file <- file.path(block_dir, filename)

      if (resume && file.exists(file)) {
        block <- readRDS(file)
        if (!identical(block$hashes, hashes) ||
            !identical(block$scenario, scenario) ||
            !identical(block$results$replicate, seq.int(start, end))) {
          stop("Checkpoint differs from current code or configuration: ", file)
        }
        rows[[length(rows) + 1L]] <- block$results
        initial <- block$initial
        next
      }

      records <- list()
      for (j in seq.int(start, end)) {
        sample <- tda_sim_generate(scenario, j, design)
        fit <- tda_test_regions(
          sample$y, design$s, design$lo, design$hi,
          work_nodes = scenario$work_nodes, initial = initial
        )
        if (fit$diagnostics$solver_ok) initial <- fit$initial
        records[[length(records) + 1L]] <-
          simulation_record(scenario, j, sample, design, fit)
      }

      result <- do.call(rbind, records)
      block <- list(
        scenario = scenario, results = result, initial = initial, hashes = hashes
      )
      temp_file <- paste0(file, ".tmp")
      saveRDS(block, temp_file)
      if (!file.rename(temp_file, file)) stop("Cannot save checkpoint: ", file)
      rows[[length(rows) + 1L]] <- result
      message(sprintf("Cell %d/276, replicates %d-%d completed", i, start, end))
    }

    results <- do.call(rbind, rows)
    write.csv(results, file.path(output_dir, sprintf("cell_%03d.csv", i)),
              row.names = FALSE)
    summaries[[length(summaries) + 1L]] <- summarize_cell(scenario, results)
    write.csv(do.call(rbind, summaries), file.path(output_dir, "summary.csv"),
              row.names = FALSE)
  }

  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  invisible(do.call(rbind, summaries))
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  cells <- if (length(arguments) >= 1L && arguments[1L] != "all") {
    as.integer(strsplit(arguments[1L], ",", fixed = TRUE)[[1L]])
  } else {
    seq_len(276L)
  }
  replicates <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 1000L
  output_dir <- if (length(arguments) >= 3L) arguments[3L] else "simulation_results"
  tda_simulation(cells, replicates, output_dir)
}
