tda_test <- function(y, s, objective = "EQ", lower_margin = -1,
                     upper_margin = 1, trial_alpha = 0.05) {
  k <- length(y)
  expand <- function(z) {
    if (!length(z) %in% c(1L, k)) {
      stop("Design inputs must have length one or length(y)")
    }
    rep_len(z, k)
  }

  objective <- expand(objective)
  objective[objective == "NI"] <- "NI_upper"
  if (anyNA(objective) ||
      any(!objective %in% c("EQ", "NI_upper", "NI_lower"))) {
    stop("objective must be EQ, NI_upper or NI_lower")
  }

  alpha <- expand(trial_alpha)
  if (any(!is.finite(alpha)) || any(alpha <= 0 | alpha >= 0.5)) {
    stop("trial_alpha must lie between 0 and 0.5")
  }

  L <- expand(lower_margin)
  U <- expand(upper_margin)
  if (!is.numeric(L) || !is.numeric(U) || anyNA(L) || anyNA(U)) {
    stop("Margins must be numeric without missing values")
  }
  if (any(!is.finite(L[objective != "NI_upper"])) ||
      any(!is.finite(U[objective != "NI_lower"]))) {
    stop("Supply finite clinical margins for the selected objective")
  }
  if (any(L[objective == "EQ"] >= U[objective == "EQ"])) {
    stop("EQ lower margins must be smaller than upper margins")
  }
  if (length(s) != k || any(!is.finite(s)) || any(s <= 0)) {
    stop("Supply matching positive standard errors")
  }

  critical <- qnorm(1 - alpha)
  lower <- ifelse(objective == "NI_upper", -Inf, L + critical * s)
  upper <- ifelse(objective == "NI_lower", Inf, U - critical * s)
  tda_test_regions(y, s, lower, upper)
}

tda_test_regions <- function(y, s, lower, upper, work_nodes = 24L,
                             initial = numeric()) {
  y <- as.numeric(y)
  s <- as.numeric(s)
  k <- length(y)

  if (k < 2L || length(s) != k || any(!is.finite(y)) ||
      any(!is.finite(s)) || any(s <= 0)) {
    stop("Supply at least two finite estimates and matching positive standard errors")
  }
  expand <- function(z) {
    if (!is.numeric(z) || !length(z) %in% c(1L, k) || anyNA(z)) {
      stop("Each endpoint must be numeric, with length one or length(y)")
    }
    rep_len(z, k)
  }
  lower <- expand(lower)
  upper <- expand(upper)
  if (!work_nodes %in% c(24L, 48L, 96L)) {
    stop("work_nodes must be 24, 48 or 96")
  }

  data <- list(y = y, s = s, lo = lower, hi = upper)
  fit <- fit_npmle(y, s)
  fit_ok <- isTRUE(fit$ok) && is.finite(fit$loglik_gap_bound) &&
    fit$loglik_gap_bound <= 1e-5

  working <- fit_working_normal(y, s)
  rule <- quadrature_rules[[as.character(work_nodes)]]
  model <- working_moments(data, working, rule)
  weights <- choose_weights(model, k, initial)
  W <- comparison_weights(model, weights)

  moments <- contrast_moments(data, fit, model, W, quadrature_rules[["24"]])
  coarse <- contrast_moments(data, fit, model, W, quadrature_rules[["12"]])
  delta <- max(abs(drop(moments$p) - drop(coarse$p)))
  refined <- FALSE
  unresolved <- FALSE

  if (delta > 1e-6 || decisions_differ(moments$p, coarse$p) ||
      any(moments$variances < -1e-14)) {
    fine <- contrast_moments(data, fit, model, W, quadrature_rules[["48"]])
    refined <- TRUE

    if (max(abs(drop(fine$p) - drop(moments$p))) > 1e-6 ||
        decisions_differ(fine$p, moments$p)) {
      finer <- contrast_moments(data, fit, model, W, quadrature_rules[["96"]])
      unresolved <- max(abs(drop(finer$p) - drop(fine$p))) > 1e-6 ||
        decisions_differ(finer$p, fine$p)
      fine <- finer
    }
    moments <- fine
  }

  solver_ok <- isTRUE(weights$ok) || isTRUE(weights$empty)
  numeric_ok <- all(is.finite(moments$V)) && all(moments$V >= 0) &&
    max(abs(moments$mean_sum)) < 1e-8 && !unresolved && solver_ok
  S <- drop(moments$S)[1L]
  V <- drop(moments$V)[1L]
  valid <- fit_ok && numeric_ok
  p <- if (valid) drop(moments$p)[1L] else 1

  if (!valid) {
    status <- "numerical_failure"
    warning(
      "Numerical checks failed; p.value is set to 1 under the simulation policy. ",
      "Inspect diagnostics.",
      call. = FALSE
    )
  } else if (weights$empty || V <= 0) {
    status <- "no_usable_cross_precision_information"
  } else {
    status <- "ok"
  }

  f <- model$f
  c <- model$c
  d <- sqrt((s[c] - s[f]) * (s[c] + s[f]))
  a <- (lower[c] - y[f]) / d
  b <- (upper[c] - y[f]) / d
  q <- ifelse(
    a > 0,
    pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE),
    pnorm(b) - pnorm(a)
  )
  A <- as.numeric(lower < upper & y > lower & y < upper)
  positive <- weights$w > 1e-10

  list(
    p.value = p,
    statistic = if (V > 0) S / sqrt(V) else NA_real_,
    S = S,
    variance = V,
    status = status,
    k = k,
    favorable = A,
    nonempty_regions = sum(lower < upper),
    informative_targets = length(unique(c[positive])),
    contributing_studies = length(unique(c(f[positive], c[positive]))),
    eligible_pairs = length(f),
    retained_pairs = length(f) - weights$excluded,
    positive_pairs = sum(positive),
    working_parameters = working,
    fit = fit,
    contrasts = data.frame(
      predictor = f, target = c, q = q, A = A[c], D = A[c] - q,
      weight = weights$w, raw_weight = weights$raw
    ),
    study_variances = moments$variances[, 1L],
    study_means = moments$means[, 1L],
    diagnostics = list(
      fit_ok = fit_ok,
      numeric_ok = numeric_ok,
      fit_gap = fit$loglik_gap_bound,
      solver_ok = solver_ok,
      empty = weights$empty,
      cap_active = weights$active,
      max_load = weights$load,
      kkt_residual = weights$residual,
      quadrature_delta = delta,
      refined = refined,
      unresolved = unresolved,
      max_variance_fraction = drop(moments$max_fraction)[1L]
    ),
    initial = weights$beta
  )
}

decisions_differ <- function(p, q) {
  any((drop(p) <= 0.10) != (drop(q) <= 0.10))
}

gauss_legendre <- function(n) {
  j <- seq_len(n - 1L)
  M <- matrix(0, n, n)
  M[cbind(j, j + 1L)] <- M[cbind(j + 1L, j)] <- j / sqrt(4 * j^2 - 1)
  e <- eigen(M, symmetric = TRUE)
  o <- order(e$values)
  list(node = e$values[o], weight = 2 * e$vectors[1L, o]^2)
}

interval_probability <- function(y, sd, lo, hi) {
  a <- (lo - y) / sd
  b <- (hi - y) / sd
  p <- ifelse(
    a > 0,
    pnorm(a, lower.tail = FALSE) - pnorm(b, lower.tail = FALSE),
    pnorm(b) - pnorm(a)
  )
  p <- pmax(0, pmin(1, p))
  p[rep_len(lo >= hi, length(p))] <- 0
  p
}

integration_grid <- function(mu, sd, lo, hi, rule) {
  ends <- (c(lo, hi) - mu) / sd
  cuts <- sort(unique(c(-9:9, ends[is.finite(ends) & ends > -9 & ends < 9])))
  half <- diff(cuts) / 2
  mid <- head(cuts, -1L) + half
  z <- as.vector(outer(rule$node, half) + rep(mid, each = length(rule$node)))
  list(y = mu + sd * z, weight = as.vector(outer(rule$weight, half)) * dnorm(z))
}

check_profile_likelihood <- function(y, s, candidate_v, tolerance = 5e-8, max_splits = 200000L) {
  ss <- s^2
  evaluate <- function(v) {
    vv <- v + ss
    mu <- sum(y / vv) / sum(1 / vv)
    rss <- sum((y - mu)^2 / vv)
    list(
      f = 0.5 * (sum(log(vv)) + rss),
      deriv = 0.5 * sum(1 / vv - (y - mu)^2 / vv^2),
      rss = rss, mu = mu
    )
  }
  original <- evaluate(candidate_v)
  zero <- evaluate(0)
  best <- min(original$f, zero$f)
  bestv <- if (zero$f < original$f) 0 else candidate_v
  interval <- function(a, b) {
    h <- (b - a) / 2
    z <- evaluate(a + h)
    zb <- evaluate(b)
    if (z$f < best) {
      best <<- z$f
      bestv <<- a + h
    }
    if (zb$f < best) {
      best <<- zb$f
      bestv <<- b
    }
    lower <- max(
      z$f - abs(z$deriv) * h - 0.25 * sum(1 / (a + ss)^2) * h^2,
      0.5 * (sum(log(a + ss)) + zb$rss)
    )
    lower <- lower - 128 * .Machine$double.eps * (1 + abs(lower) + abs(z$f))
    c(a = a, b = b, lower = lower)
  }
  right <- diff(range(y))^2
  queue <- if (right > 0) matrix(interval(0, right), nrow = 1L) else matrix(numeric(), 0L, 3L)
  splits <- 0L
  while (nrow(queue) && splits < max_splits) {
    j <- which.min(queue[, 3L])
    if (best - queue[j, 3L] <= tolerance) break
    a <- queue[j, 1L]
    b <- queue[j, 2L]
    mid <- (a + b) / 2
    if (mid == a || mid == b) break
    queue <- rbind(queue[-j, , drop = FALSE], interval(a, mid), interval(mid, b))
    splits <- splits + 1L
  }
  lower <- if (nrow(queue)) min(best, min(queue[, 3L])) else best
  list(
    mu = evaluate(bestv)$mu, tau2 = bestv,
    candidate_gap_upper = max(0, original$f - lower),
    gap = best - lower, converged = best - lower <= tolerance
  )
}

fit_working_normal <- function(y, s) {
  fn <- function(v) {
    wt <- 1 / (v + s^2)
    mu <- sum(wt * y) / sum(wt)
    0.5 * sum(log(v + s^2) + wt * (y - mu)^2)
  }
  upper <- max(1, 10 * var(y))
  fit <- optimize(fn, c(0, upper), tol = 1e-10)
  while (fit$minimum > 0.9 * upper) {
    upper <- upper * 4
    fit <- optimize(fn, c(0, upper), tol = 1e-10)
  }
  v <- if (fn(0) <= fit$objective + 1e-10) 0 else fit$minimum
  wt <- 1 / (v + s^2)
  candidate <- c(mu = sum(wt * y) / sum(wt), tau2 = v)
  audit <- check_profile_likelihood(y, s, v)
  if (!audit$converged) stop("Global working-likelihood gap check did not converge")
  if (audit$candidate_gap_upper <= 1e-7) candidate else c(mu = audit$mu, tau2 = audit$tau2)
}

update_mixture <- function(f, y, s) {
  f$weights <- f$weights / sum(f$weights)
  K <- exp(-0.5 * (outer(y, f$support, "-") / s)^2)
  f$q_scaled <- drop(K %*% f$weights)
  f$loglik <- if (all(f$q_scaled > 0)) sum(log(f$q_scaled) - log(s) - 0.5 * log(2 * pi)) else -Inf
  f
}

optimize_mixture <- function(f, y, s, maxit) {
  m <- length(f$support)
  n <- length(y)
  cached <- NULL
  calculate <- function(p) {
    if (!is.null(cached) && identical(cached$p, p)) return(cached)
    locations <- p[seq_len(m)]
    weights <- p[m + seq_len(m)]
    residual <- outer(y, locations, "-")
    K <- exp(-0.5 * (residual / s)^2)
    q <- drop(K %*% weights)
    if (any(q <= 0)) return(list(value = 1e100, gradient = rep(0, 2L * m)))
    R <- K / q
    gradient <- c(-weights * colSums(R * residual / s^2), n - colSums(R))
    cached <<- list(p = p, value = n * sum(weights) - sum(log(q)), gradient = gradient)
    cached
  }
  fit <- optim(
    c(f$support, f$weights),
    fn = function(p) calculate(p)$value,
    gr = function(p) calculate(p)$gradient,
    method = "L-BFGS-B",
    lower = c(rep(min(y), m), rep(0, m)),
    upper = c(rep(max(y), m), rep(1, m)),
    control = list(lmm = 10L, factr = 10, pgtol = 1e-9, maxit = maxit)
  )
  failures <- f$solver_failures + as.integer(fit$convergence != 0L)
  keep <- fit$par[m + seq_len(m)] > 0
  if (!any(keep)) {
    f$solver_failures <- failures
    return(f)
  }
  candidate <- f
  candidate$support <- fit$par[seq_len(m)][keep]
  candidate$weights <- fit$par[m + seq_len(m)][keep]
  candidate <- update_mixture(candidate, y, s)
  if (is.finite(candidate$loglik) && candidate$loglik >= f$loglik - 1e-9) f <- candidate
  f$solver_failures <- failures
  f
}

maximize_mixture_score <- function(f, y, s, tolerance, certify, cap = 50000L) {
  dual <- function(t) colSums(exp(-0.5 * (outer(y, t, "-") / s)^2) / f$q_scaled)
  grid <- sort(unique(c(f$support, y, seq(min(y), max(y), length.out = 65L))))
  values <- dual(grid)
  j <- which.max(values)
  best <- values[j]
  location <- grid[j]
  if (length(grid) > 2L) for (j in 2L:(length(grid) - 1L)) {
    if (values[j] >= values[j - 1L] && values[j] >= values[j + 1L]) {
      z <- optimize(function(t) -dual(t), c(grid[j - 1L], grid[j + 1L]), tol = 1e-11)
      if (-z$objective > best) {
        best <- -z$objective
        location <- z$minimum
      }
    }
  }
  if (!certify) return(list(t = location, value = best))
  bound <- function(a, b) {
    t <- pmax(a, pmin(b, y))
    envelope <- sum(exp(-0.5 * ((y - t) / s)^2) / f$q_scaled)
    h <- function(z) ifelse(abs(z) > 40, 0, abs(z^2 - 1) * exp(-0.5 * z^2))
    za <- (a - y) / s
    zb <- (b - y) / s
    hmax <- pmax(h(za), h(zb))
    for (z in c(-sqrt(3), 0, sqrt(3))) hmax <- pmax(hmax, ifelse(z >= za & z <= zb, h(z), 0))
    curvature <- sum(hmax / (s^2 * f$q_scaled))
    smooth <- max(dual(c(a, b))) + curvature * (b - a)^2 / 8
    min(envelope, smooth) * (1 + 1e-12) + 1e-12
  }
  if (length(grid) == 1L) {
    return(list(
      t = location, value = best, upper = best * (1 + 1e-12) + 1e-12,
      intervals = 0L, cap = FALSE
    ))
  }
  left <- head(grid, -1L)
  right <- tail(grid, -1L)
  queue <- cbind(left, right, vapply(seq_along(left), function(j) bound(left[j], right[j]), 0.0))
  count <- 0L
  while (nrow(queue) && count < cap) {
    j <- which.max(queue[, 3L])
    if (length(y) * log(max(best, queue[j, 3L]) / best) <= tolerance) break
    a <- queue[j, 1L]
    b <- queue[j, 2L]
    if (queue[j, 3L] <= best) {
      queue <- queue[-j, , drop = FALSE]
      next
    }
    mid <- a + (b - a) / 2
    if (mid == a || mid == b) break
    value <- dual(mid)
    if (value > best) {
      best <- value
      location <- mid
    }
    candidates <- rbind(c(a, mid, bound(a, mid)), c(mid, b, bound(mid, b)))
    queue <- rbind(queue[-j, , drop = FALSE], candidates[candidates[, 3L] >= best, , drop = FALSE])
    count <- count + 1L
  }
  list(
    t = location, value = best,
    upper = if (nrow(queue)) max(best, queue[, 3L]) else best,
    intervals = count, cap = count >= cap
  )
}

npmle_fit <- function(y, s, tolerance = 1e-5, max_outer = 80L, maxit = 1000L) {
  n <- length(y)
  quantiles <- sort(y)[floor((n - 1L) * c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 1L]
  weighted_mean <- sum(y / s^2) / sum(1 / s^2)
  atoms <- sort(unique(c(quantiles, weighted_mean)))
  f <- list(support = atoms, weights = rep(1 / length(atoms), length(atoms)), ok = FALSE,
            loglik_gap_bound = Inf, solver_failures = 0L, intervals = 0L, interval_cap = FALSE)
  f <- update_mixture(f, y, s)
  for (iter in seq_len(max_outer)) {
    f$outer <- iter
    f <- optimize_mixture(f, y, s, maxit)
    if (!is.finite(f$loglik)) break
    mx <- maximize_mixture_score(f, y, s, tolerance * 0.2, FALSE)
    if (n * log(max(1, mx$value / n)) <= tolerance * 0.8 || iter == max_outer) {
      mx <- maximize_mixture_score(f, y, s, tolerance * 0.2, TRUE)
      f$intervals <- f$intervals + mx$intervals
      f$interval_cap <- f$interval_cap || mx$cap
      f$loglik_gap_bound <- n * log(max(1, mx$upper / n))
      if (f$loglik_gap_bound <= tolerance) {
        f$ok <- TRUE
        break
      }
    }
    if (iter == max_outer) break
    z <- exp(-0.5 * ((y - mx$t) / s)^2) / f$q_scaled - 1
    a <- 0
    b <- 1 - 1e-12
    for (j in seq_len(55L)) {
      mid <- (a + b) / 2
      if (sum(z / (1 + mid * z)) > 0) a <- mid else b <- mid
    }
    alpha <- (a + b) / 2
    if (alpha < 1e-14) alpha <- 1e-10
    f$support <- c(f$support, mx$t)
    f$weights <- c(f$weights * (1 - alpha), alpha)
    f <- update_mixture(f, y, s)
  }
  f
}

fit_npmle <- function(y, s) {
  f <- npmle_fit(y, s)
  if (!f$ok) {
    g <- npmle_fit(y, s, max_outer = 160L, maxit = 2000L)
    if (g$ok || g$loglik > f$loglik) f <- g
  }
  f
}

working_moments <- function(x, par, g = quadrature_rules[["24"]]) {
  k <- length(x$s)
  edges <- which(outer(x$s, x$s, "<") & matrix(x$lo < x$hi, k, k, byrow = TRUE), arr.ind = TRUE)
  f <- edges[, 1L]
  c <- edges[, 2L]
  m <- nrow(edges)
  d <- sqrt((x$s[c] - x$s[f]) * (x$s[c] + x$s[f]))
  pi <- interval_probability(par[1L], sqrt(x$s^2 + par[2L]), x$lo, x$hi)
  Sigma <- matrix(0, m, m)
  b <- numeric(m)
  for (i in seq_len(k)) {
    incoming <- which(c == i)
    outgoing <- which(f == i)
    if (!length(incoming) && !length(outgoing)) next
    grid <- integration_grid(par[1L], sqrt(x$s[i]^2 + par[2L]), x$lo[i], x$hi[i], g)
    Q <- matrix(as.numeric(grid$y > x$lo[i] & grid$y < x$hi[i]), ncol = 1L)
    if (length(outgoing)) {
      q <- vapply(outgoing, function(e) {
        interval_probability(grid$y, d[e], x$lo[c[e]], x$hi[c[e]])
      }, numeric(length(grid$y)))
      Q <- cbind(Q, q)
    }
    mean <- c(pi[i], pi[c[outgoing]])
    C <- crossprod(Q * sqrt(grid$weight)) - tcrossprod(mean)
    C[1L, 1L] <- pi[i] * (1 - pi[i])
    if (length(incoming)) {
      Sigma[incoming, incoming] <-
        Sigma[incoming, incoming, drop = FALSE] + C[1L, 1L]
    }
    if (length(outgoing)) {
      b[outgoing] <- pi[c[outgoing]] * (1 - pi[c[outgoing]]) - C[1L, -1L]
      Sigma[outgoing, outgoing] <-
        Sigma[outgoing, outgoing, drop = FALSE] + C[-1L, -1L, drop = FALSE]
      if (length(incoming)) {
        cross <- matrix(C[-1L, 1L], length(outgoing), length(incoming))
        Sigma[outgoing, incoming] <- Sigma[outgoing, incoming, drop = FALSE] - cross
        Sigma[incoming, outgoing] <- Sigma[incoming, outgoing, drop = FALSE] - t(cross)
      }
    }
  }
  list(Sigma = Sigma, b = b, f = f, c = c, d = d, pi = pi)
}

solve_nonnegative_qp <- function(Q, a, initial = numeric(), tolerance = 1e-8, maxit = 120L) {
  n <- length(a)
  beta <- if (length(initial) == n) pmax(0, initial) else numeric(n)
  mask <- beta > 0
  if (!any(mask)) mask[] <- TRUE
  threshold <- tolerance * max(1, max(abs(a)))
  best <- n + 1L
  budget <- 3L
  for (iter in seq_len(maxit)) {
    beta[] <- 0
    if (any(mask)) {
      R <- chol(Q[mask, mask, drop = FALSE])
      beta[mask] <- backsolve(R, backsolve(R, a[mask], transpose = TRUE))
    }
    gradient <- drop(Q %*% beta) - a
    bad <- which((mask & beta < -1e-12) | (!mask & gradient < -threshold))
    if (!length(bad)) {
      beta <- pmax(0, beta)
      gradient <- drop(Q %*% beta) - a
      residual <- max(ifelse(beta > 0, abs(gradient), pmax(0, -gradient)))
      if (residual <= threshold) break
      bad <- which(beta == 0 & gradient < -threshold)
      if (!length(bad)) break
    }
    if (length(bad) < best) {
      best <- length(bad)
      budget <- 3L
    } else budget <- budget - 1L
    if (budget >= 0L) mask[bad] <- !mask[bad] else mask[tail(bad, 1L)] <- !mask[tail(bad, 1L)]
  }
  beta <- pmax(0, beta)
  gradient <- drop(Q %*% beta) - a
  residual <- max(ifelse(beta > 0, abs(gradient), pmax(0, -gradient)))
  list(
    beta = beta, residual = residual, ok = residual <= threshold,
    active = FALSE, iterations = iter
  )
}

solve_load_constrained_qp <- function(Q, a, P, cap = 10, tol = 1e-8,
                                      maxit = 20000L, start = numeric()) {
  n <- length(a)
  k <- nrow(P)
  beta <- if (length(start) == n) pmax(0, start) else numeric(n)
  maxload <- max(P %*% beta)
  if (maxload > cap) beta <- beta * cap / maxload * (1 - 32 * .Machine$double.eps)
  free <- which(beta > 0)
  active <- integer()
  multiplier <- numeric(k)
  threshold <- tol * max(1, max(abs(a)))
  for (iter in seq_len(maxit)) {
    gradient <- drop(Q %*% beta) - a
    direction <- numeric(n)
    multiplier[] <- 0
    if (length(free)) {
      Qf <- Q[free, free, drop = FALSE]
      if (length(active)) {
        Pf <- P[active, free, drop = FALSE]
        K <- rbind(cbind(Qf, t(Pf)), cbind(Pf, matrix(0, length(active), length(active))))
        rhs <- c(-gradient[free], rep(0, length(active)))
        solution <- tryCatch(solve(K, rhs), error = function(e) qr.solve(K, rhs, tol = 1e-12))
        direction[free] <- solution[seq_along(free)]
        multiplier[active] <- tail(solution, length(active))
      } else direction[free] <- solve(Qf, -gradient[free])
    }
    if (max(abs(direction)) > 1e-11 * (1 + max(abs(beta)))) {
      step <- 1
      blocking_variable <- blocking_cap <- NA_integer_
      decreasing <- free[direction[free] < -1e-14]
      if (length(decreasing)) {
        limits <- -beta[decreasing] / direction[decreasing]
        j <- which.min(limits)
        if (limits[j] < step) {
          step <- max(0, limits[j])
          blocking_variable <- decreasing[j]
        }
      }
      loads <- drop(P %*% beta)
      changes <- drop(P %*% direction)
      increasing <- setdiff(which(changes > 1e-14), active)
      if (length(increasing)) {
        limits <- (cap - loads[increasing]) / changes[increasing]
        j <- which.min(limits)
        if (limits[j] < step) {
          step <- max(0, limits[j])
          blocking_cap <- increasing[j]
          blocking_variable <- NA_integer_
        }
      }
      beta <- beta + step * direction
      beta[abs(beta) < 1e-14] <- 0
      if (!is.na(blocking_variable)) {
        beta[blocking_variable] <- 0
        free <- setdiff(free, blocking_variable)
      }
      if (!is.na(blocking_cap)) active <- c(active, blocking_cap)
      next
    }
    reduced <- gradient + drop(crossprod(P, multiplier))
    zero <- setdiff(seq_len(n), free)
    bad_zero <- if (length(zero)) zero[which.min(reduced[zero])] else NA_integer_
    bad_cap <- if (length(active)) active[which.min(multiplier[active])] else NA_integer_
    vz <- if (is.na(bad_zero)) 0 else min(0, reduced[bad_zero])
    vc <- if (is.na(bad_cap)) 0 else min(0, multiplier[bad_cap])
    if (min(vz, vc) >= -threshold) break
    if (vc < vz) active <- setdiff(active, bad_cap) else free <- c(free, bad_zero)
  }
  beta <- pmax(0, beta)
  loads <- drop(P %*% beta)
  if (max(loads) > cap) beta <- beta * cap / max(loads) * (1 - 32 * .Machine$double.eps)
  loads <- drop(P %*% beta)
  gradient <- drop(Q %*% beta) - a
  reduced <- gradient + drop(crossprod(P, multiplier))
  residual <- max(
    ifelse(beta > 1e-10, abs(reduced), pmax(0, -reduced)),
    pmax(0, loads - cap), pmax(0, -multiplier),
    abs(multiplier * (cap - loads))
  )
  list(
    beta = beta, residual = residual,
    ok = is.finite(residual) && residual < 1e-7 * max(1, max(abs(a))),
    active = any(multiplier > 1e-10), iterations = iter
  )
}

study_loads <- function(model, raw, k) {
  v <- numeric(k)
  for (e in seq_along(raw)) {
    v[model$f[e]] <- v[model$f[e]] + raw[e]
    v[model$c[e]] <- v[model$c[e]] + raw[e]
  }
  v
}

choose_weights <- function(model, k, initial = numeric()) {
  sd <- sqrt(pmax(0, diag(model$Sigma)))
  use <- which(sd^2 >= 1e-12)
  m <- length(sd)
  if (!length(use)) {
    return(list(
      ok = FALSE, empty = TRUE, w = numeric(m), raw = numeric(m),
      beta = numeric(m), load = 0, active = FALSE,
      residual = NA_real_, excluded = m
    ))
  }
  Q <- model$Sigma[use, use, drop = FALSE] / outer(sd[use], sd[use]) + diag(0.01, length(use))
  a <- model$b[use] / sd[use]
  P <- matrix(0, k, length(use))
  P[cbind(model$f[use], seq_along(use))] <- 1 / sd[use]
  P[cbind(model$c[use], seq_along(use))] <- 1 / sd[use]
  fit <- solve_nonnegative_qp(Q, a, if (length(initial) == m) initial[use] else numeric())
  if (!fit$ok || max(P %*% fit$beta) > 10) {
    fit <- solve_load_constrained_qp(Q, a, P, start = fit$beta)
  }
  raw <- beta <- numeric(m)
  beta[use] <- fit$beta
  raw[use] <- fit$beta / sd[use]
  list(
    ok = fit$ok && sum(raw) > 0, empty = FALSE,
    w = if (sum(raw) > 0) raw / sum(raw) else raw,
    raw = raw, beta = beta, load = max(study_loads(model, raw, k)),
    active = fit$active, residual = fit$residual, excluded = m - length(use)
  )
}

comparison_weights <- function(model, opt) {
  m <- length(model$f)
  if (!m) return(matrix(numeric(), 0L, 3L))
  count <- tabulate(model$c, nbins = max(model$c))
  target <- 1 / count[model$c]
  cbind(tda = opt$w, equal = rep(1 / m, m), target = target / sum(target))
}

contrast_moments <- function(x, fit, model, W, g = quadrature_rules[["24"]]) {
  k <- length(x$y)
  J <- ncol(W)
  f <- model$f
  c <- model$c
  d <- model$d
  pi <- vapply(seq_len(k), function(i) {
    sum(fit$weights * interval_probability(
      fit$support, x$s[i], x$lo[i], x$hi[i]
    ))
  }, 0.0)
  means <- variances <- observed <- matrix(0, k, J)
  for (i in seq_len(k)) {
    incoming <- which(c == i)
    outgoing <- which(f == i)
    wi <- if (length(incoming)) colSums(W[incoming, , drop = FALSE]) else numeric(J)
    gi <- function(y) {
      value <- outer(as.numeric(y > x$lo[i] & y < x$hi[i]), wi)
      if (length(outgoing)) {
        q <- vapply(outgoing, function(e) {
          interval_probability(y, d[e], x$lo[c[e]], x$hi[c[e]])
        }, numeric(length(y)))
        q <- matrix(q, nrow = length(y))
        value <- value - q %*% W[outgoing, , drop = FALSE]
      }
      value
    }
    observed[i, ] <- gi(x$y[i])
    predicted <- if (length(outgoing)) {
      drop(crossprod(W[outgoing, , drop = FALSE], pi[c[outgoing]]))
    } else {
      numeric(J)
    }
    means[i, ] <- wi * pi[i] - predicted
    if (!length(outgoing)) {
      variances[i, ] <- wi^2 * pi[i] * (1 - pi[i])
    } else for (h in seq_along(fit$support)) {
      if (fit$weights[h] <= 0) next
      grid <- integration_grid(fit$support[h], x$s[i], x$lo[i], x$hi[i], g)
      centered <- sweep(gi(grid$y), 2L, means[i, ], "-")
      variances[i, ] <- variances[i, ] + fit$weights[h] * colSums(centered^2 * grid$weight)
    }
  }
  S <- colSums(observed)
  V <- colSums(variances)
  p <- ifelse(V > 0, pnorm(S / sqrt(V), lower.tail = FALSE), 1)
  fraction <- ifelse(V > 0, apply(variances, 2L, max) / V, 0)
  list(
    S = S, V = V, p = p, mean_sum = colSums(means),
    max_fraction = fraction, means = means, variances = variances
  )
}

egger_p <- function(y, s, design) {
  if (max(s) - min(s) < 1e-12) return(NA_real_)
  X <- cbind(1, 1 / s)
  fit <- lm.fit(X, y / s)
  df <- length(y) - 2L
  if (fit$rank < 2L || df <= 0L) return(NA_real_)
  covariance <- chol2inv(qr.R(fit$qr))
  se <- sqrt(sum(fit$residuals^2) / df * covariance[1L, 1L])
  statistic <- fit$coefficients[1L] / se
  if (design == "NI") pt(statistic, df) else 2 * pt(-abs(statistic), df)
}

quadrature_rules <- setNames(
  lapply(c(12L, 24L, 48L, 96L), gauss_legendre),
  c("12", "24", "48", "96")
)
