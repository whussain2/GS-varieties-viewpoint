# ==================================================================
# Fig. 3C - favorable-allele frequencies under recurrent selection
# JXB Viewpoint, Waseem Hussain (PlantNura)
# 20 inbred parents, 200 additive QTL on 5 chromosomes, 200 DH lines per cycle
# from 20 crosses, 20 parents selected per cycle on a predictor of accuracy 0.6.
# Records the frequency of the favorable allele at every QTL in every cycle.
# Writes Fig3C_locus_frequencies.csv
# ==================================================================
set.seed(3)
n_chr <- 5; n_loci <- 240; n_tot <- n_chr * n_loci; chr <- rep(1:n_chr, each = n_loci)
n_qtl <- 200; n_par <- 20; n_cand <- 200; n_cycles <- 12; r_acc <- 0.6

p <- runif(n_tot, 0.15, 0.85)
par <- matrix(0, n_par, n_tot); for (j in 1:n_tot) par[, j] <- 2 * rbinom(n_par, 1, p[j])
qtl <- sort(sample(1:n_tot, n_qtl)); eff <- rgamma(n_qtl, 0.6, 1) * sample(c(-1, 1), n_qtl, TRUE)
fav <- ifelse(eff > 0, 2, 0)
freq <- matrix(NA, n_cycles + 1, n_qtl); best <- mean_fav <- numeric(n_cycles + 1)
for (cyc in 0:n_cycles) {
  pop <- matrix(0, n_cand, n_tot); row <- 0
  for (cr in 1:20) {
    pair <- sample(1:n_par, 2)
    for (k in 1:10) {
      row <- row + 1
      for (h in 1:n_chr) {
        loci <- which(chr == h); from <- rep(sample(1:2, 1), n_loci)
        for (x in sort(sample(1:n_loci, rpois(1, 1)))) from[x:n_loci] <- 3 - from[x:n_loci]
        pop[row, loci] <- ifelse(from == 1, par[pair[1], loci], par[pair[2], loci])
      }
    }
  }
  isfav <- sweep(pop[, qtl], 2, fav, "==")
  freq[cyc + 1, ] <- colMeans(isfav); nf <- rowSums(isfav); best[cyc + 1] <- max(nf); mean_fav[cyc + 1] <- mean(nf)
  g <- as.vector(pop[, qtl] %*% eff)
  pred <- r_acc * scale(g)[, 1] + sqrt(1 - r_acc^2) * rnorm(n_cand)
  par <- pop[order(pred, decreasing = TRUE)[1:n_par], ]
}
out <- data.frame(cycle = rep(0:n_cycles, n_qtl), locus = rep(1:n_qtl, each = n_cycles + 1),
                  abs_effect = rep(abs(eff), each = n_cycles + 1), freq = as.vector(freq))
write.csv(out, "Fig3C_locus_frequencies.csv", row.names = FALSE)
write.csv(data.frame(cycle = 0:n_cycles, mean_freq = rowMeans(freq), best_line = best, mean_line = mean_fav), "Fig3C_summary.csv", row.names = FALSE)
cat("mean freq:", round(rowMeans(freq), 3), "\nbest line (of 200):", best, "\n")
cat("fixed favorable at cycle 12:", sum(freq[13, ] == 1), " lost:", sum(freq[13, ] == 0), "\n")
