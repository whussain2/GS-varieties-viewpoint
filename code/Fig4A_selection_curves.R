# ==================================================================
# Fig. 4A data: selection on a criterion of accuracy r   (JXB Viewpoint, Waseem Hussain, PlantNura)
#
# Panel B : what selecting the top 10% keeps, for ANY selection criterion
#           of accuracy r (genomic or phenotypic) - simulation, base R
# Panel C : why variety counts lag: population mean under GS and under
#           phenotypic selection (from the 20-year simulation, Fig. 3)
#           and the genetic value of varieties reaching release, which
#           trails the population by the product-development lag
#
# Inputs : Fig3_multicycle_summary.csv (written by Fig3_multicycle_GS_vs_PS.R)
# Packages: ggplot2
# ==================================================================

library(ggplot2)
set.seed(2026)

teal <- "#0F766E"; orange <- "#C2410C"; purple <- "#6D28D9"; slate <- "#334155"; mute <- "#94A3B8"
base_family <- "Inter"
theme_box <- theme_minimal(base_size = 8.5, base_family = base_family) +
  theme(plot.title = element_text(face = "bold", size = 9.5, colour = slate),
        plot.subtitle = element_text(size = 7.6, colour = "#475569", margin = margin(b = 4)),
        axis.title = element_text(size = 8, colour = slate),
        axis.text = element_text(size = 7, colour = "#475569"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#E2E8F0", linewidth = 0.3),
        legend.position = "top", legend.justification = "left",
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        legend.key.height = unit(7, "pt"), legend.key.width = unit(14, "pt"), legend.spacing.y = unit(0, "pt"),
        plot.margin = margin(6, 8, 6, 6))

# ------------------------------------------------------------------
# Panel B: select the top 10% of 1000 candidates on a criterion with accuracy r
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

pB <- ggplot(dB, aes(r, p, colour = what)) +
  annotate("rect", xmin = 0.66, xmax = 0.80, ymin = -Inf, ymax = Inf, fill = "#F1F5F9") +
  geom_vline(xintercept = c(0.66, 0.80), colour = c(teal, slate), linetype = "22", linewidth = 0.4) +
  annotate("text", x = 0.65, y = 0.42, label = "GEBV,\nuntested", hjust = 1, size = 2.4, colour = teal,
           family = base_family, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 0.81, y = 0.62, label = "field\nmeans", hjust = 0, size = 2.4, colour = slate,
           family = base_family, fontface = "bold", lineheight = 0.9) +
  geom_hline(yintercept = c(0.5, 0.1), colour = mute, linetype = "dotted", linewidth = 0.35) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  scale_colour_manual(values = c(slate, teal, orange)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Accuracy of the selection criterion, r", y = "Probability / share",
       title = NULL, subtitle = NULL) +
  theme_box

