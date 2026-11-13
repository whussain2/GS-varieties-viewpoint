# ------------------------------------------------------------------
# Repeat the Figure 2 simulation 100 times
#
# Figure 2 shows ONE simulated breeding programme. One run can be lucky
# or unlucky, so here we repeat the whole thing (new parents, new crosses,
# new field trial, new GBLUP fit) 100 times and average the results.
# The numbers quoted in the Comment come from this script.
#
# The steps inside the loop are the same as in Fig2_GBLUP_simulation.R,
# written out again so each script runs on its own.
# Plain base R + rrBLUP (only to estimate the two variances by REML).
# ------------------------------------------------------------------

library(rrBLUP)
set.seed(2026)

n_rep  <- 100
result <- data.frame()

for (rep_i in 1:n_rep) {

  # ---- 1. parents and loci ----
  n_par <- 20; n_chr <- 5; n_loci <- 240; n_tot <- n_chr * n_loci
  chr <- rep(1:n_chr, each = n_loci)
  p   <- runif(n_tot, 0.1, 0.9)
  parents <- matrix(0, n_par, n_tot)
  for (j in 1:n_tot) parents[, j] <- 2 * rbinom(n_par, 1, p[j])
  qtl_pos    <- sort(sample(1:n_tot, 200))
  marker_pos <- setdiff(1:n_tot, qtl_pos)
  qtl_effect <- rnorm(200, 0, 1)

  # ---- 2. RILs: 30 training crosses + 20 candidate crosses, 10 RILs each ----
  n_cross <- 50; ril_per_cross <- 10; n_lines <- 500
  geno <- matrix(0, n_lines, n_tot)
  cross <- rep(1:n_cross, each = ril_per_cross)
  pairs_all <- t(combn(n_par, 2)); pick <- sample(1:nrow(pairs_all), n_cross)   # 50 different pairs
  row <- 0
  for (c in 1:n_cross) {
    pair <- pairs_all[pick[c], ]
    for (k in 1:ril_per_cross) {
      row <- row + 1
      for (h in 1:n_chr) {
        loci_h <- which(chr == h)
        xo   <- sort(sample(1:n_loci, rpois(1, 2)))
        from <- rep(sample(1:2, 1), n_loci)
        for (x in xo) from[x:n_loci] <- 3 - from[x:n_loci]
        geno[row, loci_h] <- ifelse(from == 1, parents[pair[1], loci_h], parents[pair[2], loci_h])
      }
    }
  }
  train <- which(cross <= 30)
  cand  <- which(cross > 30)
  g_true <- as.vector(geno[, qtl_pos] %*% qtl_effect)
  g_true <- (g_true - mean(g_true[train])) / sd(g_true[train])

  # ---- 3. field trial: 2 env x 2 reps, plot h2 = 0.30 ----
  se2_plot <- (1 - 0.30) / 0.30
  plots <- matrix(0, length(train), 4)
  for (b in 1:4) plots[, b] <- 50 + c(-2, -2, 3, 3)[b] + c(0.5, -0.5, 0.3, -0.3)[b] +
                                g_true[train] + rnorm(length(train), 0, sqrt(se2_plot))
  plots <- sweep(plots, 2, colMeans(plots))        # remove block effects
  y <- rowMeans(plots)                             # adjusted line means

  # ---- 4. G matrix ----
  M <- geno[, marker_pos]
  freq <- colMeans(M) / 2
  keep <- freq > 0.02 & freq < 0.98
  M <- M[, keep]; freq <- freq[keep]
  W <- sweep(M, 2, 2 * freq)
  G <- W %*% t(W) / (2 * sum(freq * (1 - freq))) + diag(0.01, n_lines)

  # ---- 5. GBLUP by hand (variances by REML from rrBLUP) ----
  vc  <- mixed.solve(y, K = G[train, train])
  sg2 <- vc$Vu; se2 <- vc$Ve; lambda <- se2 / sg2
  n <- length(y)
  X <- matrix(1, n, 1)
  Z <- matrix(0, n, n_lines); Z[cbind(1:n, train)] <- 1
  LHS <- rbind(cbind(t(X) %*% X, t(X) %*% Z),
               cbind(t(Z) %*% X, t(Z) %*% Z + lambda * solve(G)))
  RHS <- rbind(t(X) %*% y, t(Z) %*% y)
  C   <- solve(LHS)
  gebv <- (C %*% RHS)[-1]
  pev  <- diag(C)[-1] * se2
  rel  <- 1 - pev / (sg2 * diag(G))

  # ---- results for this replicate ----
  wf <- c()                                   # accuracy within each candidate cross (family)
  for (cc in 31:50) {
    ii <- which(cross == cc)
    wf <- c(wf, cor(gebv[ii], g_true[ii]))
  }
  pheno_dev <- y - mean(y)
  sel  <- cand[order(gebv[cand],  decreasing = TRUE)[1:20]]   # top 10% on GEBV
  top  <- cand[order(g_true[cand], decreasing = TRUE)[1:20]]  # truly top 10%
  best <- cand[which.max(g_true[cand])]
  sel_gain <- mean(g_true[sel]) - mean(g_true[cand])          # true gain from selecting on GEBV
  max_gain <- mean(g_true[top]) - mean(g_true[cand])          # gain with perfect knowledge

  result <- rbind(result, data.frame(
    rep = rep_i,
    h2_line_mean      = sg2 * mean(diag(G[train, train])) / (sg2 * mean(diag(G[train, train])) + se2),
    shrinkage_slope   = unname(coef(lm(gebv[train] ~ pheno_dev))[2]),
    rel_training      = mean(rel[train]),
    rel_candidates    = mean(rel[cand]),
    acc_training      = cor(gebv[train], g_true[train]),
    acc_candidates    = cor(gebv[cand],  g_true[cand]),
    acc_within_family = mean(wf),
    pheno_acc_training = cor(y, g_true[train]),
    top20_overlap     = length(intersect(sel, top)),
    best_in_top10pct  = best %in% sel,
    top_gebv_is_best  = cand[which.max(gebv[cand])] == best,
    above_avg_in_sel  = mean(g_true[sel] > mean(g_true[cand])),
    gain_vs_perfect   = sel_gain / max_gain))
  cat("replicate", rep_i, "done\n")
}

write.csv(result, "Box1_100_replicates.csv", row.names = FALSE)

# 95% binomial confidence intervals for the two yes/no outcomes
print(binom.test(sum(result$best_in_top10pct), n_rep)$conf.int)
print(binom.test(sum(result$top_gebv_is_best), n_rep)$conf.int)

# ---- summary: mean, sd and 5-95% range over the 100 replicates ----
num <- result[, -1]
num$best_in_top10pct <- as.numeric(num$best_in_top10pct)
num$top_gebv_is_best <- as.numeric(num$top_gebv_is_best)
summary_rep <- data.frame(quantity = names(num),
                          mean = round(colMeans(num), 3),
                          sd   = round(apply(num, 2, sd), 3),
                          q05  = round(apply(num, 2, quantile, 0.05), 3),
                          q95  = round(apply(num, 2, quantile, 0.95), 3))
print(summary_rep, row.names = FALSE)
write.csv(summary_rep, "Box1_100_replicates_summary.csv", row.names = FALSE)
