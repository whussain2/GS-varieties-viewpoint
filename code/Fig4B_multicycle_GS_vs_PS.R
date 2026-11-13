# ------------------------------------------------------------------
# Fig. 4 (JXB Viewpoint): genetic gain over 20 years
#   PS          : phenotypic selection, one cycle every 4 years
#   GS-retrain  : genomic selection, one cycle per year, training set
#                 updated every year with new field data (2-year lag)
#   GS-static   : genomic selection, one cycle per year, model trained
#                 once at year 0 and never updated
#
# The three schemes use the SAME field budget: 200 plots per year.
#   PS  : 200 lines x 4 plots, once every 4 years  = 200 plots/year
#   GS  : 50 lines  x 4 plots, every year           = 200 plots/year
#   (GS also pays for genotyping 200 lines per year.)
#
# Inbred crop with doubled-haploid (DH) lines; purely additive trait.
# Questions addressed:
#   - does rapid cycling still win once genetic variance and accuracy
#     are allowed to erode?
#   - what happens if the training set is not updated?
#   - how does GS compare with phenotypic selection on the same budget?
#
# Plain base R. rrBLUP (mixed.solve) is used only to fit the marker model.
# ------------------------------------------------------------------

library(rrBLUP)
set.seed(2026)

n_rep   <- 30          # replicates of the whole 20-year programme
n_years <- 20
n_par   <- 20          # parents selected each cycle
n_cross <- 20          # crosses per cycle
n_dh    <- 10          # DH lines per cross -> 200 candidates
n_chr   <- 5; n_loci <- 240; n_tot <- n_chr * n_loci
chr     <- rep(1:n_chr, each = n_loci)
h2_plot <- 0.30        # plot-level heritability in the base population
n_plots <- 4           # plots per phenotyped line (2 env x 2 reps)

gain_PS <- matrix(NA, n_rep, n_years + 1)   # mean true value of the current population
gain_GR <- matrix(NA, n_rep, n_years + 1)
gain_GS <- matrix(NA, n_rep, n_years + 1)
sd_PS   <- matrix(NA, n_rep, n_years + 1)   # genetic SD of the current population
sd_GR   <- matrix(NA, n_rep, n_years + 1)
sd_GS   <- matrix(NA, n_rep, n_years + 1)
acc_GR  <- matrix(NA, n_rep, n_years)       # accuracy of GEBVs among candidates
acc_GS  <- matrix(NA, n_rep, n_years)

for (rep_i in 1:n_rep) {

  # ---------- founders, QTL and markers ----------
  p <- runif(n_tot, 0.1, 0.9)
  founders <- matrix(0, n_par, n_tot)
  for (j in 1:n_tot) founders[, j] <- 2 * rbinom(n_par, 1, p[j])
  qtl  <- sort(sample(1:n_tot, 200))
  mrk  <- setdiff(1:n_tot, qtl)
  eff  <- rnorm(200)

  # base population: 200 DH lines from 20 crosses among the founders
  par_geno <- founders
  base <- matrix(0, n_cross * n_dh, n_tot); row <- 0
  for (c in 1:n_cross) {
    pair <- sample(1:nrow(par_geno), 2)
    for (k in 1:n_dh) {
      row <- row + 1
      for (h in 1:n_chr) {
        loci <- which(chr == h)
        from <- rep(sample(1:2, 1), n_loci)
        for (x in sort(sample(1:n_loci, rpois(1, 1)))) from[x:n_loci] <- 3 - from[x:n_loci]  # one meiosis
        base[row, loci] <- ifelse(from == 1, par_geno[pair[1], loci], par_geno[pair[2], loci])
      }
    }
  }
  g_base  <- as.vector(base[, qtl] %*% eff)
  mu0     <- mean(g_base); sd0 <- sd(g_base)          # everything is scaled to the base population
  se_plot <- sqrt(sd0^2 * (1 - h2_plot) / h2_plot)    # residual SD per plot, fixed for all years
  y_base  <- g_base + rnorm(length(g_base), 0, se_plot / sqrt(n_plots))   # line means

  # ================= PS: one cycle every 4 years =================
  pop <- base; y <- y_base
  gain_PS[rep_i, 1] <- 0; sd_PS[rep_i, 1] <- 1
  for (yr in 1:n_years) {
    if (yr %% 4 == 0) {                                # a new cycle finishes every 4 years
      parents <- pop[order(y, decreasing = TRUE)[1:n_par], ]
      new <- matrix(0, n_cross * n_dh, n_tot); row <- 0
      for (c in 1:n_cross) {
        pair <- sample(1:n_par, 2)
        for (k in 1:n_dh) {
          row <- row + 1
          for (h in 1:n_chr) {
            loci <- which(chr == h)
            from <- rep(sample(1:2, 1), n_loci)
            for (x in sort(sample(1:n_loci, rpois(1, 1)))) from[x:n_loci] <- 3 - from[x:n_loci]
            new[row, loci] <- ifelse(from == 1, parents[pair[1], loci], parents[pair[2], loci])
          }
        }
      }
      pop <- new
      g   <- as.vector(pop[, qtl] %*% eff)
      y   <- g + rnorm(length(g), 0, se_plot / sqrt(n_plots))   # all 200 lines phenotyped
    }
    g <- as.vector(pop[, qtl] %*% eff)
    gain_PS[rep_i, yr + 1] <- (mean(g) - mu0) / sd0
    sd_PS[rep_i, yr + 1]   <- sd(g) / sd0
  }

  # ================= GS: one cycle per year (two versions) ========
  for (scheme in c("retrain", "static")) {
    pop   <- base
    train_M <- base[, mrk]; train_y <- y_base          # year-0 training set
    pending <- list()                                  # field data that arrive later
    gain <- c(0, rep(NA, n_years)); gsd <- c(1, rep(NA, n_years)); acc <- rep(NA, n_years)

    # fit the year-0 model (used every year by the static scheme)
    fit0 <- mixed.solve(train_y, Z = train_M)
    for (yr in 1:n_years) {
      # 1. predict the current population and choose parents
      if (scheme == "retrain") {
        fit <- mixed.solve(train_y, Z = train_M)
      } else {
        fit <- fit0
      }
      gebv <- as.vector(pop[, mrk] %*% fit$u)
      g    <- as.vector(pop[, qtl] %*% eff)
      acc[yr] <- cor(gebv, g)
      parents <- pop[order(gebv, decreasing = TRUE)[1:n_par], ]

      # 2. phenotype 50 random lines of the current population (results arrive in 2 years)
      ph  <- sample(1:nrow(pop), 50)
      pending[[yr]] <- list(M = pop[ph, mrk], y = g[ph] + rnorm(50, 0, se_plot / sqrt(n_plots)))
      if (yr > 2 && scheme == "retrain") {             # add data from 2 years ago, keep last 4 cohorts + base
        train_M <- base[, mrk]; train_y <- y_base
        for (k in max(1, yr - 5):(yr - 2)) {
          train_M <- rbind(train_M, pending[[k]]$M); train_y <- c(train_y, pending[[k]]$y)
        }
      }

      # 3. make the next population (DH lines, one year per cycle with speed breeding)
      new <- matrix(0, n_cross * n_dh, n_tot); row <- 0
      for (c in 1:n_cross) {
        pair <- sample(1:n_par, 2)
        for (k in 1:n_dh) {
          row <- row + 1
          for (h in 1:n_chr) {
            loci <- which(chr == h)
            from <- rep(sample(1:2, 1), n_loci)
            for (x in sort(sample(1:n_loci, rpois(1, 1)))) from[x:n_loci] <- 3 - from[x:n_loci]
            new[row, loci] <- ifelse(from == 1, parents[pair[1], loci], parents[pair[2], loci])
          }
        }
      }
      pop <- new
      g   <- as.vector(pop[, qtl] %*% eff)
      gain[yr + 1] <- (mean(g) - mu0) / sd0
      gsd[yr + 1]  <- sd(g) / sd0
    }
    if (scheme == "retrain") { gain_GR[rep_i, ] <- gain; sd_GR[rep_i, ] <- gsd; acc_GR[rep_i, ] <- acc }
    if (scheme == "static")  { gain_GS[rep_i, ] <- gain; sd_GS[rep_i, ] <- gsd; acc_GS[rep_i, ] <- acc }
  }
  cat("replicate", rep_i, "done\n")
}

# ---------------- summaries ----------------
years <- 0:n_years
out <- data.frame(year = years,
  PS_gain = colMeans(gain_PS), PS_se = apply(gain_PS, 2, sd) / sqrt(n_rep),
  GSretrain_gain = colMeans(gain_GR), GSretrain_se = apply(gain_GR, 2, sd) / sqrt(n_rep),
  GSstatic_gain = colMeans(gain_GS), GSstatic_se = apply(gain_GS, 2, sd) / sqrt(n_rep),
  PS_genSD = colMeans(sd_PS), GSretrain_genSD = colMeans(sd_GR), GSstatic_genSD = colMeans(sd_GS),
  GSretrain_acc = c(NA, colMeans(acc_GR, na.rm = TRUE)), GSstatic_acc = c(NA, colMeans(acc_GS, na.rm = TRUE)))
# (accuracy is NA in a replicate once all genetic variance is lost)
print(round(out[c(1, 5, 9, 11, 13, 17, 21), ], 3))
write.csv(round(out, 4), "Fig4B_multicycle_summary.csv", row.names = FALSE)

# The figure is drawn by Fig4_plot.R
