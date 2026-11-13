# ==================================================================
# Fig. 4A data: selection on a criterion of accuracy r
# JXB Viewpoint, Waseem Hussain (PlantNura)
#
# What selecting the top 10% of 1000 candidates keeps, for ANY selection
# criterion of accuracy r (genomic or phenotypic). Plain base R.
# Writes Fig4A_selection_curves.csv (plotted by Fig4_plot.R)
# ==================================================================

set.seed(2026)

# ------------------------------------------------------------------
# select the top 10% of 1000 candidates on a criterion with accuracy r
# ------------------------------------------------------------------
N <- 1000; n_keep <- 100; n_rep <- 2000
acc <- seq(0.05, 0.95, by = 0.05)
best_in <- top1 <- above <- numeric(length(acc))
for (k in 1:length(acc)) {
  r <- acc[k]; q1 <- q2 <- q3 <- 0
  for (rep in 1:n_rep) {
    g    <- rnorm(N)
    pred <- r * g + sqrt(1 - r^2) * rnorm(N)
    kept <- order(pred, decreasing = TRUE)[1:n_keep]
    best <- which.max(g)
    q1 <- q1 + (best %in% kept)
    q2 <- q2 + (which.max(pred) == best)
    q3 <- q3 + mean(g[kept] > 0)
  }
  best_in[k] <- q1 / n_rep; top1[k] <- q2 / n_rep; above[k] <- q3 / n_rep
}
dB <- rbind(data.frame(r = acc, p = above,   what = "Selected lines truly above average"),
            data.frame(r = acc, p = best_in, what = "Truly best line is selected"),
            data.frame(r = acc, p = top1,    what = "Top-ranked line is the truly best"))
dB$what <- factor(dB$what, levels = unique(dB$what))
write.csv(dB, "Fig4A_selection_curves.csv", row.names = FALSE)
