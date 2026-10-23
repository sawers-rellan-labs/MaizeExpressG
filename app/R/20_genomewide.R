# Genome-wide operations, vectorised over genes.

# Design matrix for the screen models.
screen_design <- function(x, model, spline_df = 3) switch(model,
  linear    = cbind(1, x),
  quadratic = cbind(1, x, x^2),
  cubic     = cbind(1, x, x^2, x^3),
  spline    = cbind(1, splines::ns(x, df = spline_df)))

# Fit the same linear model to every gene at once.
#   Y: genes x libraries (already transformed), s: matching sample table,
#   x: positions for the columns of Y.
screen_genes <- function(Y, s, x, model = "linear", spline_df = 3) {
  n <- ncol(Y)
  X <- screen_design(x, model, spline_df)
  qx <- qr(X); p <- qx$rank
  if (n - p < 1) stop("Not enough libraries/sections in the selected range for this model.")
  B <- t(qr.coef(qx, t(Y)))
  Fit <- t(qr.fitted(qx, t(Y)))
  rss <- rowSums((Y - Fit)^2)
  tss <- rowSums((Y - rowMeans(Y))^2)
  fstat <- ((tss - rss) / (p - 1)) / (rss / (n - p))
  out <- data.frame(id = rownames(Y), R2 = 1 - rss / tss,
                    F = fstat, p = pf(fstat, p - 1, n - p, lower.tail = FALSE),
                    stringsAsFactors = FALSE)
  if (model == "linear") {
    se <- sqrt(rss / (n - p) * solve(crossprod(X))[2, 2])
    out$slope <- B[, 2]; out$slope_se <- se
  }
  # Section means -> peak position and end-vs-start change
  secs <- sort(unique(s$section))
  A <- sapply(secs, function(k) (s$section == k) / sum(s$section == k))
  M <- Y %*% A
  xs <- tapply(x, s$section, `[`, 1)[as.character(secs)]
  out$peak_x <- xs[max.col(M, ties.method = "first")]
  out$end_minus_start <- M[, ncol(M)] - M[, 1]
  # One-way ANOVA across sections and lack-of-fit (need replicate libraries)
  k <- length(secs)
  ss_pe <- rowSums((Y - M %*% t(A > 0))^2)
  if (n - k > 0) {
    f_sec <- ((tss - ss_pe) / (k - 1)) / (ss_pe / (n - k))
    out$section_p <- pf(f_sec, k - 1, n - k, lower.tail = FALSE)
    if (k - p > 0) {
      f_lof <- ((rss - ss_pe) / (k - p)) / (ss_pe / (n - k))
      out$lackfit_p <- pf(f_lof, k - p, n - k, lower.tail = FALSE)
    }
  }
  out$q <- p.adjust(out$p, "BH")
  out
}

# Correlate a query profile (one value per section) with every gene's section means.
correlate_profiles <- function(M, q, method = "pearson") {
  ok <- is.finite(q)
  M <- M[, ok, drop = FALSE]; q <- q[ok]
  if (method == "spearman") {
    M <- t(apply(M, 1, rank)); q <- rank(q)
  }
  Mc <- M - rowMeans(M); qc <- q - mean(q)
  r <- as.vector(Mc %*% qc) / sqrt(rowSums(Mc^2) * sum(qc^2))
  k <- length(q)
  tt <- r * sqrt((k - 2) / pmax(1 - r^2, 1e-12))
  p <- 2 * pt(-abs(tt), k - 2)
  data.frame(id = rownames(M), r = r, p = p, q = p.adjust(p, "BH"), stringsAsFactors = FALSE)
}

# Template shapes for pattern matching (k sections).
template_profile <- function(kind, k, at) switch(kind,
  increasing = seq_len(k),
  decreasing = rev(seq_len(k)),
  peak       = exp(-((seq_len(k) - at)^2) / 4),
  dip        = -exp(-((seq_len(k) - at)^2) / 4),
  step_up    = as.numeric(seq_len(k) >= at),
  step_down  = as.numeric(seq_len(k) < at))

# Correlation between maize and rice orthologs on the UDM unified gradient.
ortholog_pairs <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    o <- LG$orthologs
    mz <- LG$udm$maize; rc <- LG$udm$rice
    o <- o[o$maize %in% rownames(mz) & o$rice %in% rownames(rc), ]
    a <- log2(mz[o$maize, , drop = FALSE] + 1); b <- log2(rc[o$rice, , drop = FALSE] + 1)
    ac <- a - rowMeans(a); bc <- b - rowMeans(b)
    r <- rowSums(ac * bc) / sqrt(rowSums(ac^2) * rowSums(bc^2))
    cache <<- data.frame(maize = o$maize, maize_symbol = GENES$symbol[match(o$maize, GENES$id)],
                         rice = o$rice, r_UDM = r,
                         maize_mean_log2 = rowMeans(a), rice_mean_log2 = rowMeans(b),
                         maize_peak = max.col(a, "first"), rice_peak = max.col(b, "first"),
                         maize_description = GENES$description[match(o$maize, GENES$id)],
                         rice_description = GENES$description[match(o$rice, GENES$id)],
                         stringsAsFactors = FALSE)
    cache
  }
})

# PCA of libraries (GEO) for QC.
library_pca <- function(sp, min_reads = 0, ntop = 2000) {
  s <- dataset_samples(sp, "geo", min_reads)
  Y <- log2(LG$expr[[sp]][, s$sample] + 1)
  Y <- Y[rowMeans(Y) > 1, ]
  v <- apply(Y, 1, var)
  Y <- Y[order(v, decreasing = TRUE)[seq_len(min(ntop, nrow(Y)))], ]
  pc <- prcomp(t(Y), center = TRUE, scale. = FALSE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  list(scores = cbind(s, PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3]), var = ve)
}

# Correlation of each library with the mean of its section (log2, expressed genes).
library_qc <- function(sp) {
  s <- dataset_samples(sp, "geo", 0)
  Y <- log2(LG$expr[[sp]][, s$sample] + 1)
  Y <- Y[rowMeans(LG$expr[[sp]][, s$sample]) > 5, ]
  r <- vapply(seq_len(nrow(s)), function(i) {
    others <- s$sample[s$section == s$section[i] & s$sample != s$sample[i]]
    if (!length(others)) return(NA_real_)
    cor(Y[, s$sample[i]], rowMeans(Y[, others, drop = FALSE]))
  }, 0)
  s$r_vs_other_reps <- r
  s
}
