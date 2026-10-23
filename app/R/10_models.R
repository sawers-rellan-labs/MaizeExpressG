# Model fitting for one gene's profile along the gradient.
# `d` has columns x (position), y (expression on the chosen scale) and section.

MODELS <- c("Linear" = "linear", "Quadratic" = "quadratic", "Cubic" = "cubic",
            "Natural spline" = "spline", "Segmented (1 breakpoint)" = "segmented",
            "Logistic (4-parameter sigmoid)" = "logistic")

model_formula_text <- function(model, spline_df = 3) switch(model,
  linear    = "y ~ x",
  quadratic = "y ~ x + I(x^2)",
  cubic     = "y ~ x + I(x^2) + I(x^3)",
  spline    = sprintf("y ~ splines::ns(x, df = %d)", spline_df),
  segmented = "y ~ x + pmax(x - c, 0)   (c = breakpoint, grid search)",
  logistic  = "y ~ SSfpl(x, A, B, xmid, scal)   (nls)")

fit_segmented <- function(d) {
  ux <- sort(unique(d$x))
  if (length(ux) < 4) return(NULL)
  grid <- seq(ux[2], ux[length(ux) - 1], length.out = 200)
  rss <- vapply(grid, function(cc) sum(resid(lm(y ~ x + pmax(x - cc, 0), data = d))^2), 0)
  cbest <- grid[which.min(rss)]
  fit <- lm(y ~ x + pmax(x - cbest, 0), data = d)
  # profile-likelihood interval for the breakpoint (F-based, 1 extra parameter)
  n <- nrow(d); p <- 4
  thr <- min(rss) * (1 + qf(0.95, 1, max(n - p, 1)) / max(n - p, 1))
  inside <- grid[rss <= thr]
  list(fit = fit, c = cbest, c_lo = min(inside), c_hi = max(inside), k = 4)
}

# 4-parameter logistic y = A + (B - A) / (1 + exp((xmid - x) / scal)), A = value at
# the base end, B = value at the tip end. Bounded fit (plateaus within half a data
# range of the observed values, midpoint inside the sampled range) from several starts.
fit_logistic <- function(d) {
  m <- tapply(d$y, d$x, mean); xs <- as.numeric(names(m))
  ry <- diff(range(d$y)); rx <- diff(range(xs))
  if (ry == 0) return(NULL)
  lower <- c(A = min(d$y) - ry / 2, B = min(d$y) - ry / 2, xmid = min(xs), scal = rx / 200)
  upper <- c(A = max(d$y) + ry / 2, B = max(d$y) + ry / 2, xmid = max(xs), scal = rx)
  best <- NULL
  for (q in c(0.25, 0.5, 0.75)) {
    st <- list(A = m[[1]], B = m[[length(m)]], xmid = quantile(xs, q, names = FALSE), scal = rx / 10)
    f <- tryCatch(nls(y ~ A + (B - A) / (1 + exp((xmid - x) / scal)), data = d, start = st,
                      algorithm = "port", lower = lower, upper = upper,
                      control = nls.control(maxiter = 500)), error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  if (!is.null(best)) attr(best, "bounds") <- list(lower = lower, upper = upper)
  best
}

fit_one <- function(d, model, spline_df = 3) {
  out <- tryCatch({
    switch(model,
      linear    = list(fit = lm(y ~ x, data = d), k = 2),
      quadratic = list(fit = lm(y ~ x + I(x^2), data = d), k = 3),
      cubic     = list(fit = lm(y ~ x + I(x^2) + I(x^3), data = d), k = 4),
      spline    = list(fit = lm(y ~ splines::ns(x, df = spline_df), data = d), k = spline_df + 1),
      segmented = fit_segmented(d),
      logistic  = list(fit = fit_logistic(d), k = 4))
  }, error = function(e) NULL)
  if (is.null(out) || is.null(out$fit)) return(NULL)
  out$k <- out$k + 1                                     # + residual variance, for AIC/BIC
  out$model <- model
  out
}

# Goodness of fit and tests for one fitted model.
model_stats <- function(f, d) {
  n <- nrow(d); y <- d$y
  rss <- sum(resid(f$fit)^2); tss <- sum((y - mean(y))^2)
  p <- f$k - 1                                            # mean-function parameters
  r2 <- 1 - rss / tss
  adj <- 1 - (1 - r2) * (n - 1) / max(n - p, 1)
  ll <- -n / 2 * (log(2 * pi) + log(rss / n) + 1)
  fstat <- ((tss - rss) / (p - 1)) / (rss / (n - p))
  p_model <- if (n > p && p > 1) pf(fstat, p - 1, n - p, lower.tail = FALSE) else NA
  # Lack-of-fit against the section-means (one-way ANOVA) model; needs replicates.
  k_sec <- length(unique(d$section))
  ss_pe <- sum(tapply(y, d$section, function(v) sum((v - mean(v))^2)))
  df_pe <- n - k_sec; df_lof <- k_sec - p
  p_lof <- if (df_pe > 0 && df_lof > 0 && ss_pe > 0)
    pf(((rss - ss_pe) / df_lof) / (ss_pe / df_pe), df_lof, df_pe, lower.tail = FALSE) else NA
  data.frame(Model = names(MODELS)[MODELS == f$model], n = n, Parameters = as.integer(p),
             RSS = rss, Sigma = sqrt(rss / max(n - p, 1)), R2 = r2, Adj_R2 = adj,
             AIC = -2 * ll + 2 * f$k, BIC = -2 * ll + log(n) * f$k,
             F_p = p_model, LackOfFit_p = p_lof, stringsAsFactors = FALSE)
}

# Prediction grid with confidence band (lm models; segmented is conditional on c).
predict_curve <- function(f, d, level = 0.95, npts = 200) {
  grid <- data.frame(x = seq(min(d$x), max(d$x), length.out = npts))
  if (inherits(f$fit, "nls")) {
    grid$fit <- predict(f$fit, newdata = grid); grid$lwr <- NA; grid$upr <- NA
  } else {
    pr <- predict(f$fit, newdata = grid, interval = "confidence", level = level)
    grid$fit <- pr[, "fit"]; grid$lwr <- pr[, "lwr"]; grid$upr <- pr[, "upr"]
  }
  grid$model <- names(MODELS)[MODELS == f$model]
  grid
}

coef_table <- function(f, level = 0.95) {
  cf <- summary(f$fit)$coefficients
  ci <- tryCatch(suppressMessages(
    if (inherits(f$fit, "nls")) confint.default(f$fit, level = level) else confint(f$fit, level = level)),
    error = function(e) NULL)
  out <- data.frame(Term = rownames(cf), Estimate = cf[, 1], SE = cf[, 2], t = cf[, 3], p = cf[, 4],
                    stringsAsFactors = FALSE)
  if (!is.null(ci)) { out$CI_low <- ci[, 1]; out$CI_high <- ci[, 2] }
  out$Term <- sub("pmax\\(x - cbest, 0\\)", "slope change after breakpoint", out$Term)
  if (f$model == "segmented")
    out <- rbind(out, data.frame(Term = "breakpoint c (x units)", Estimate = f$c, SE = NA, t = NA, p = NA,
                                 CI_low = f$c_lo, CI_high = f$c_hi))
  rownames(out) <- NULL
  out
}

# Plain-language quantities derived from the focus model.
derived_quantities <- function(f, d, scale, xunit) {
  cf <- coef(f$fit); lines <- character()
  per <- paste0("per ", xunit)
  if (f$model == "linear") {
    b <- cf[["x"]]
    lines <- c(lines, sprintf("Slope: %.4g %s %s.", b, if (scale == "log2") "log2 units" else "FPKM", per))
    if (scale == "log2") lines <- c(lines, sprintf(
      "= %.3g-fold change %s; expression %s every %.3g %s.",
      2^b, per, if (b > 0) "doubles" else "halves", 1 / abs(b), xunit))
  }
  if (f$model == "quadratic") {
    v <- -cf[[2]] / (2 * cf[[3]])
    kind <- if (cf[[3]] < 0) "maximum" else "minimum"
    inr <- v >= min(d$x) && v <= max(d$x)
    lines <- c(lines, sprintf("Vertex (%s) at x = %.3g %s%s.", kind, v, xunit,
                              if (inr) "" else " - outside the fitted range"))
  }
  if (f$model %in% c("cubic", "spline")) {
    g <- predict_curve(f, d)
    lines <- c(lines, sprintf("Fitted curve peaks at x = %.3g and is lowest at x = %.3g (%s).",
                              g$x[which.max(g$fit)], g$x[which.min(g$fit)], xunit))
  }
  if (f$model == "segmented") {
    lines <- c(lines,
      sprintf("Breakpoint at x = %.3g %s (approx. 95%% profile CI %.3g-%.3g).", f$c, xunit, f$c_lo, f$c_hi),
      sprintf("Slope before: %.4g; after: %.4g (%s).", cf[[2]], cf[[2]] + cf[[3]], per))
  }
  if (f$model == "logistic") {
    lines <- c(lines,
      sprintf("Base-end plateau A = %.4g, tip-end plateau B = %.4g (%s along the gradient).",
              cf[["A"]], cf[["B"]], if (cf[["B"]] > cf[["A"]]) "rising" else "falling"),
      sprintf("Half-way transition at xmid = %.3g %s; 10-90%% transition width = %.3g %s.",
              cf[["xmid"]], xunit, 2 * log(9) * cf[["scal"]], xunit))
    b <- attr(f$fit, "bounds")
    at <- names(cf)[abs(cf - b$lower[names(cf)]) < 1e-6 | abs(cf - b$upper[names(cf)]) < 1e-6]
    if (length(at)) lines <- c(lines, sprintf(
      "Parameter(s) %s hit the fitting bounds - a plateau or the transition may lie outside the sampled range.",
      paste(at, collapse = ", ")))
  }
  lines
}

# Nested F-tests between polynomial models and against the section-means model.
nested_tests <- function(d) {
  m1 <- lm(y ~ x, d); m2 <- lm(y ~ x + I(x^2), d); m3 <- lm(y ~ x + I(x^2) + I(x^3), d)
  out <- list()
  a <- anova(m1, m2); out[[1]] <- c("Quadratic vs linear", a$F[2], a$`Pr(>F)`[2])
  a <- anova(m2, m3); out[[2]] <- c("Cubic vs quadratic", a$F[2], a$`Pr(>F)`[2])
  if (any(duplicated(d$section))) {
    m0 <- lm(y ~ 1, d); mf <- lm(y ~ factor(section), d)
    a <- anova(m0, mf); out[[3]] <- c("Any section effect (one-way ANOVA)", a$F[2], a$`Pr(>F)`[2])
  }
  out <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  names(out) <- c("Test", "F", "p")
  out$F <- as.numeric(out$F); out$p <- as.numeric(out$p)
  out
}

# Do several genes share the same slope? (ANCOVA, linear in x)
slope_comparison <- function(d) {
  d$gene <- factor(d$gene)
  if (nlevels(d$gene) < 2) return(NULL)
  full <- lm(y ~ gene * x, d); red <- lm(y ~ gene + x, d)
  a <- anova(red, full)
  per <- lm(y ~ 0 + gene + gene:x, d)
  cf <- summary(per)$coefficients
  sl <- cf[grepl(":x$", rownames(cf)), , drop = FALSE]
  list(p = a$`Pr(>F)`[2], F = a$F[2],
       slopes = data.frame(Gene = sub("^gene(.*):x$", "\\1", rownames(sl)), Slope = sl[, 1],
                           SE = sl[, 2], p = sl[, 4], row.names = NULL))
}
