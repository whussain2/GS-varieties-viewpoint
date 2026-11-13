# ==================================================================
# Fig. 1 (Box 1) - one simulated breeding programme, start to finish
# JXB Viewpoint, Waseem Hussain (PlantNura)
#
# 1  8 bi-parental crosses -> F1 -> F2 ... F5 by single-seed descent (240 F5 lines)
# 2  genotype all lines at 600 SNPs (0 / 1 / 2)
# 3  build G (VanRaden) and choose a training set of 80 lines by CDmean
# 4  grow the training set at 3 locations x 2 replicates; record yield (t/ha)
# 5  fit RR-BLUP (marker effects) and GBLUP (breeding values); BayesB for comparison
# 6  rank all 240 lines on GEBV, with PEV and reliability
#
# Writes Fig1_data.json (values shown in the figure) and CSV files.
# Packages: rrBLUP, BGLR, jsonlite
# ==================================================================
library(rrBLUP); library(BGLR); library(jsonlite)
set.seed(11)

n_cross <- 8; n_per <- 30; N <- n_cross * n_per; n_chr <- 6; n_snp <- 600; per_chr <- n_snp / n_chr
chr <- rep(1:n_chr, each = per_chr)
pos <- unlist(lapply(1:n_chr, function(h) sort(runif(per_chr, 0, 1.5))))      # Morgans

# ---- parents: 16 inbreds simulated from 4 founder groups; the crossing plan below uses eight of them
p0  <- runif(n_snp, 0.1, 0.9)
grp <- rep(1:4, each = 4)
pg  <- sapply(1:4, function(k) pmin(pmax(p0 + rnorm(n_snp, 0, 0.3), 0.02), 0.98))
par <- t(sapply(1:16, function(i) 2 * rbinom(n_snp, 1, pg[, grp[i]])))
pairs <- cbind(seq(1, 16, 2), seq(2, 16, 2))[sample(8), ]
pairs[, 2] <- c(5, 9, 13, 1, 15, 3, 7, 11)[1:8]   # crosses between groups: eight (odd-numbered) parents, each used in two crosses

# ---- F5 lines by single-seed descent: mosaic of the two parents + residual heterozygosity (1/16)
X <- matrix(0, N, n_snp); cross <- rep(1:n_cross, each = n_per)
for (i in 1:N) {
  a <- par[pairs[cross[i], 1], ]; b <- par[pairs[cross[i], 2], ]
  for (h in 1:n_chr) {
    loci <- which(chr == h); xo <- runif(rpois(1, 1.5 * 1.8), 0, 1.5)       # RIL map expansion ~1.8
    src <- (sapply(pos[loci], function(p) sum(xo < p)) + sample(0:1, 1)) %% 2
    X[i, loci] <- ifelse(src == 0, a[loci], b[loci])
  }
  het <- which(a != b & runif(n_snp) < 1 / 16); X[i, het] <- 1
}
lines <- sprintf("G%03d", 1:N); rownames(X) <- lines; colnames(X) <- sprintf("SNP%03d", 1:n_snp)
poly <- apply(X, 2, var) > 0; chr_kept <- chr[poly]
X <- X[, poly]; n_snp <- ncol(X); colnames(X) <- sprintf("SNP%03d", 1:n_snp)   # 538 polymorphic SNPs, renumbered

# ---- true genetic values: 150 QTL among the SNPs, yield SD 0.45 t/ha
qtl <- sample(1:n_snp, 150); a_eff <- rnorm(150)
g <- as.vector(scale(X[, qtl] %*% a_eff)) * 0.45

# ---- G matrix (VanRaden 2008)
p <- colMeans(X) / 2; W <- sweep(X, 2, 2 * p); G <- tcrossprod(W) / (2 * sum(p * (1 - p)))

# ---- training set by CDmean (Rincent et al. 2012), greedy forward selection, lambda = 1
lam <- 1; n_train <- 80
cdmean <- function(T) { K <- G[, T] %*% solve(G[T, T] + lam * diag(length(T))); mean(rowSums(K * G[, T]) / diag(G)) }
T <- integer(0)
while (length(T) < n_train) {
  cand <- setdiff(1:N, T)
  if (length(T) < 2) cand <- sample(cand, 60)
  sc <- sapply(cand, function(c) cdmean(c(T, c))); T <- c(T, cand[which.max(sc)])
}
train <- sort(T); cand_set <- setdiff(1:N, train)
pc <- prcomp(W); pv <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2), 1)

# ---- multi-location trial: 3 locations x 2 reps (RCBD)
locs <- c("Location 1", "Location 2", "Location 3"); loc_eff <- c(0, -0.6, 0.5)
fld <- expand.grid(line = train, loc = 1:3, rep = 1:2)
ge  <- matrix(rnorm(N * 3, 0, 0.3), N, 3)
fld$yield <- 5 + loc_eff[fld$loc] + c(0.1, -0.1)[fld$rep] + g[fld$line] + ge[cbind(fld$line, fld$loc)] + rnorm(nrow(fld), 0, 0.9)
fld$yield <- round(fld$yield, 2)
fit <- lm(yield ~ factor(loc) / factor(rep) + factor(line), data = fld)
adj <- tapply(fld$yield - (fitted(fit) - predict(fit, newdata = transform(fld, loc = 1, rep = 1)) ), fld$line, mean)
# simpler, transparent adjusted mean: line mean over the balanced trial (locations and reps are balanced)
adj <- tapply(fld$yield, fld$line, mean)
y <- as.vector(adj); names(y) <- lines[train]

# ---- RR-BLUP: y = 1 mu + W beta + e
rr <- mixed.solve(y, Z = W[train, ]); beta <- as.vector(rr$u); gebv_rr <- as.vector(W %*% beta)
# ---- GBLUP: y = 1 mu + Z g + e, g ~ N(0, G sg2), all 240 lines in G
Z <- matrix(0, length(train), N); Z[cbind(1:length(train), train)] <- 1
gb <- mixed.solve(y, Z = Z, K = G, SE = TRUE)
gebv <- as.vector(gb$u); pev <- as.vector(gb$u.SE)^2; sg2 <- gb$Vu; se2 <- gb$Ve
rel <- pmax(0, 1 - pev / (diag(G) * sg2)); mu <- as.vector(gb$beta)
h2_mean <- sg2 / (sg2 + se2)
# ---- BayesB for comparison
bb <- BGLR(y = y, ETA = list(list(X = W[train, ], model = "BayesB")), nIter = 6000, burnIn = 1000, verbose = FALSE, saveAt = tempfile())
beta_bb <- bb$ETA[[1]]$b; gebv_bb <- as.vector(W %*% beta_bb)

cat("cor RR vs GBLUP:", round(cor(gebv, gebv_rr), 4), " cor BayesB vs GBLUP:", round(cor(gebv, gebv_bb), 3), "\n")
cat("accuracy, untested:", round(cor(gebv[cand_set], g[cand_set]), 3), " training:", round(cor(gebv[train], g[train]), 3), "\n")
cat("sg2", round(sg2, 3), "se2", round(se2, 3), "h2(line mean)", round(h2_mean, 2), "mu", round(mu, 2), "\n")
cat("CDmean selected:", round(cdmean(train), 3), " random 80:", round(mean(replicate(20, cdmean(sample(1:N, 80)))), 3), "\n")

# ---- values to display: 6 training lines from different crosses
# six training lines from six different crosses, spanning the range of GEBVs
ord <- train[order(gebv[train], decreasing = TRUE)]; show_tr <- integer(0)
for (q in c(0.03, 0.97, 0.25, 0.7, 0.45, 0.88)) {
  pool <- ord[!(cross[ord] %in% cross[show_tr])]
  show_tr <- c(show_tr, pool[which.min(abs(match(pool, ord) / length(ord) - q))])
}
# ten SNPs that vary among those lines, with the largest RR-BLUP effects (keep at least two with heterozygotes)
vary <- which(apply(X[show_tr, ], 2, function(v) length(unique(v)) > 1))
het  <- vary[apply(X[show_tr, vary], 2, function(v) any(v == 1))]
show_snp <- head(het[order(-abs(beta[het]))], 2)
show_snp <- sort(c(show_snp, head(setdiff(vary[order(-abs(beta[vary]))], show_snp), 8)))
ex_i <- train[which.max(y)]
ranked <- order(gebv, decreasing = TRUE)
top <- ranked[1:24]
out <- list(
  n = N, n_train = length(train), n_snp = n_snp, sg2 = sg2, se2 = se2, h2 = h2_mean, mu = mu, pv = pv,
  acc_untested = cor(gebv[cand_set], g[cand_set]), cor_rr = cor(gebv, gebv_rr), cor_bb = cor(gebv, gebv_bb),
  pairs = pairs, lines_show = lines[show_tr], cross_show = cross[show_tr],
  X_show = X[show_tr, show_snp], snp_show = colnames(X)[show_snp],
  G_show = round(G[show_tr, show_tr], 2), y_show = round(y[lines[show_tr]], 2),
  beta_show = round(beta[show_snp], 4), beta_bb_show = round(beta_bb[show_snp], 4),
  gebv_show = round(gebv[show_tr], 2), gebv_rr_show = round(gebv_rr[show_tr], 2),
  field = fld[fld$line %in% show_tr, ], locs = locs,
  pca = data.frame(pc1 = pc$x[, 1], pc2 = pc$x[, 2], cross = cross, train = (1:N) %in% train),
  beta_all = data.frame(snp = 1:n_snp, chr = chr_kept, beta = beta),
  rank = data.frame(rank = 1:N, line = lines[ranked], gebv = gebv[ranked], se = sqrt(pev[ranked]), rel = rel[ranked],
                    tested = ranked %in% train, true_g = g[ranked]),
  shrink = data.frame(line = lines[train], dev = y - mean(y), gebv = gebv[train]),
  ex = list(line = lines[ex_i], ybar = y[lines[ex_i]], gebv = gebv[ex_i], mean_y = mean(y), rel = rel[ex_i])
)
out$field$line <- lines[out$field$line]; out$field$loc <- locs[out$field$loc]
write_json(out, "Fig1_data.json", digits = 5, auto_unbox = TRUE)
write.csv(data.frame(line = lines, cross = cross, train = (1:N) %in% train, adj_mean = ifelse((1:N) %in% train, y[lines], NA),
                     gebv = gebv, pev = pev, reliability = rel, gebv_rrblup = gebv_rr, gebv_bayesB = gebv_bb, true_g = g),
          "Fig1_lines.csv", row.names = FALSE)
write.csv(X, "Fig1_markers_0_1_2.csv"); write.csv(round(G, 4), "Fig1_G_matrix.csv"); write.csv(fld, "Fig1_field_trial.csv", row.names = FALSE)
write.csv(data.frame(snp = colnames(X), rrblup = beta, bayesB = beta_bb), "Fig1_marker_effects.csv", row.names = FALSE)
