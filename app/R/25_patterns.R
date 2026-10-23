# Opposite base-to-tip patterns, shared by the "Base vs tip" tab and scripts/05_gradient_patterns.R.
#   base_high_declining: highest section in the basal third, tip within `tip_fold` of the lowest section
#   base_low_rising:     lowest section in the basal third, tip within `tip_fold` of the highest section
# plus >= `min_fold` change between base and tip (first / last `end_frac` of sections, log2 means)
# and a linear trend along the leaf (FDR < `fdr`, same model as the Gradient screen).

PATTERNS <- c("High at base, declining" = "base_high_declining",
              "Low at base, rising"     = "base_low_rising")

pattern_ids <- function(res) res[[intersect(c("v3_id", "gene_id"), names(res))[1]]]

find_patterns <- function(sp, dataset = "geo", min_reads = 0, min_fpkm = 5, min_fold = 4,
                          fdr = 0.01, tip_fold = 2, end_frac = 0.2) {
  s <- dataset_samples(sp, dataset, min_reads)
  lin <- section_mean_matrix(sp, dataset, "linear", min_reads)
  lin <- lin[apply(lin, 1, max) >= min_fpkm, , drop = FALSE]
  M <- section_mean_matrix(sp, dataset, "log2", min_reads, ids = rownames(lin))
  lin <- lin[rownames(M), , drop = FALSE]
  k <- ncol(M); pos <- tapply(s$pos_cm, s$section, `[`, 1)
  n_end <- max(2, round(k * end_frac))
  base <- seq_len(n_end); tip <- (k - n_end + 1):k

  fit <- screen_genes(log2(dataset_matrix(sp, dataset)[rownames(M), s$sample, drop = FALSE] + 1),
                      s, s$pos_cm, "linear")
  fit <- fit[match(rownames(M), fit$id), ]

  base_l <- rowMeans(M[, base, drop = FALSE]); tip_l <- rowMeans(M[, tip, drop = FALSE])
  lfc <- tip_l - base_l
  hi <- apply(M, 1, max); lo <- apply(M, 1, min)
  max_sec <- max.col(M, "first"); min_sec <- max.col(-M, "first")
  basal <- seq_len(ceiling(k / 3)); tol <- log2(tip_fold); need <- log2(min_fold)
  sig <- !is.na(fit$q) & fit$q < fdr
  pattern <- ifelse(lfc <= -need & max_sec %in% basal & tip_l - lo <= tol & fit$slope < 0 & sig,
                    "base_high_declining",
             ifelse(lfc >= need & min_sec %in% basal & tip_l - hi >= -tol & fit$slope > 0 & sig,
                    "base_low_rising", NA_character_))
  rk <- t(apply(M, 1, rank)); rkc <- rk - rowMeans(rk); xc <- seq_len(k) - (k + 1) / 2
  rho <- round(as.vector(rkc %*% xc) / sqrt(rowSums(rkc^2) * sum(xc^2)), 3)   # Spearman; rounded so
                                                                              # exact 0.8 ties classify stably
  # where along the leaf half of the base-to-tip change is reached
  mid <- (base_l + tip_l) / 2
  half <- vapply(seq_len(nrow(M)), function(i) {
    j <- if (lfc[i] > 0) which(M[i, ] >= mid[i])[1] else which(M[i, ] <= mid[i])[1]
    pos[[j]]
  }, 0)

  g <- GENES[match(rownames(M), GENES$id), ]
  out <- data.frame(v5_id = g$v5, v3_id = rownames(M), v3_to_v5 = g$v5_mapping, symbol = g$symbol,
                    description = g$description, pattern = pattern,
                    shape = ifelse(abs(rho) >= 0.8, "gradual", "early switch, then plateau"),
                    base_FPKM = signif(rowMeans(lin[, base, drop = FALSE]), 4),
                    tip_FPKM = signif(rowMeans(lin[, tip, drop = FALSE]), 4),
                    log2FC_tip_vs_base = round(lfc, 3), slope_log2_per_cm = signif(fit$slope, 4),
                    q_linear = signif(fit$q, 3), spearman_rho = rho, half_change_cm = half,
                    stringsAsFactors = FALSE, row.names = NULL)
  if (!any(g$species == "maize")) {                          # rice: no v3/v5 columns
    out$v5_id <- out$v3_to_v5 <- NULL
    names(out)[names(out) == "v3_id"] <- "gene_id"
  }
  means <- as.data.frame(signif(lin, 4), row.names = NULL)
  names(means) <- sprintf("FPKM_%.1fcm", pos)
  out <- cbind(out, means)
  out <- out[order(-abs(out$log2FC_tip_vs_base)), ]
  rownames(out) <- NULL
  attr(out, "z") <- t(scale(t(M)))                          # z-scored log2 section means
  attr(out, "pos") <- pos
  attr(out, "settings") <- list(species = sp, dataset = dataset, min_reads = min_reads, min_fpkm = min_fpkm,
                                min_fold = min_fold, fdr = fdr, tip_fold = tip_fold, end_frac = end_frac,
                                n_base_tip = n_end, n_sections = k)
  out
}

# All genes of each pattern as z-scored profiles with median and 10-90% band.
pattern_figure <- function(res, xlab = "Distance from leaf base (cm)", ligule_x = NULL, title = NULL,
                           subtitle = NULL, caption = NULL) {
  z <- attr(res, "z"); pos <- attr(res, "pos")
  keep <- pattern_ids(res)[!is.na(res$pattern)]
  if (!length(keep)) return(NULL)
  z <- z[keep, , drop = FALSE]
  long <- data.frame(id = rep(rownames(z), ncol(z)), x = rep(as.numeric(pos), each = nrow(z)),
                     z = as.vector(z))
  long$pattern <- factor(res$pattern[match(long$id, pattern_ids(res))], PATTERNS)
  n <- table(factor(res$pattern, PATTERNS))
  lab <- setNames(sprintf("%s\n(n = %s)", names(PATTERNS), format(as.integer(n), big.mark = ",")), PATTERNS)
  summ <- do.call(rbind, lapply(split(long, list(long$pattern, long$x), drop = TRUE), function(d)
    data.frame(pattern = d$pattern[1], x = d$x[1], med = median(d$z),
               lo = quantile(d$z, 0.1), hi = quantile(d$z, 0.9))))
  summ$pattern <- factor(summ$pattern, PATTERNS)
  ggplot(long, aes(x, z)) +
    geom_line(aes(group = id), colour = "grey80", linewidth = 0.15, alpha = 0.4) +
    geom_ribbon(data = summ, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "#2E7D32", alpha = 0.2) +
    geom_line(data = summ, aes(x, med), colour = "#1B5E20", linewidth = 1.2) +
    { if (!is.null(ligule_x)) geom_vline(xintercept = ligule_x, linetype = "dotted", colour = "grey40") } +
    facet_wrap(~pattern, labeller = labeller(pattern = lab), drop = FALSE) +
    labs(x = xlab, y = "z-score of log2(FPKM + 1)", title = title, subtitle = subtitle, caption = caption) +
    theme_minimal(base_size = 12)
}
