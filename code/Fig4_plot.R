# ==================================================================
# Fig. 4 - what genomic selection delivers
# (A) selecting on an imperfect predictor keeps a good set, not the best line
#     (data from Fig4A_selection_curves.R -> Fig4A_selection_curves.csv)
# (B) genetic gain over 20 years: GS with early recycling vs phenotypic selection
#     (data from Fig4B_multicycle_GS_vs_PS.R -> Fig4B_multicycle_summary.csv)
# Packages: ggplot2, patchwork
# ==================================================================
library(ggplot2); library(patchwork)
fam <- "Inter"
th <- theme_minimal(base_size = 8, base_family = fam) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E2E8F0", linewidth = 0.3),
        axis.title = element_text(size = 8, colour = "#334155"), axis.text = element_text(size = 7, colour = "#64748B"),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(size = 7.3), legend.key.width = unit(16, "pt"), legend.key.height = unit(8, "pt"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(14, 8, 4, 4),
        plot.title = element_text(face = "bold", size = 8.6, colour = "#0F172A"))

# ---- A: selection curves
a <- read.csv("Fig4A_selection_curves.csv")
a$what <- factor(a$what, levels = c("Selected lines truly above average", "Truly best line is selected",
                                    "Top-ranked line is the truly best"))
pA <- ggplot(a, aes(r, p, colour = what)) +
  annotate("rect", xmin = 0.66, xmax = 0.80, ymin = -Inf, ymax = Inf, fill = "#F1F5F9") +
  geom_vline(xintercept = 0.66, colour = "#0F766E", linetype = "22", linewidth = 0.4) +
  geom_vline(xintercept = 0.80, colour = "#334155", linetype = "22", linewidth = 0.4) +
  annotate("text", x = 0.645, y = 0.40, label = "GEBV,\nuntested", hjust = 1, size = 2.5, colour = "#0F766E",
           fontface = "bold", family = fam, lineheight = 0.9) +
  annotate("text", x = 0.815, y = 0.62, label = "field\nmeans", hjust = 0, size = 2.5, colour = "#334155",
           fontface = "bold", family = fam, lineheight = 0.9) +
  geom_hline(yintercept = c(0.5, 0.1), colour = "#94A3B8", linetype = "dotted", linewidth = 0.35) +
  geom_line(linewidth = 0.8) + geom_point(size = 1) +
  scale_colour_manual(values = c("#334155", "#0F766E", "#C2410C")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Accuracy of the selection criterion, r", y = "Probability / share",
       title = "A prediction keeps a good set, not the best line") + th

# ---- B: gain, two schemes
d <- read.csv("Fig4B_multicycle_summary.csv")
b <- rbind(data.frame(year = d$year, y = d$GSretrain_gain, se = d$GSretrain_se, s = "GS, 1-year cycle (early recycling)"),
           data.frame(year = d$year, y = d$PS_gain, se = d$PS_se, s = "Phenotypic selection, 4-year cycle"))
b$s <- factor(b$s, levels = unique(b$s))
pB <- ggplot(b, aes(year, y, colour = s, fill = s)) +
  geom_ribbon(data = b[b$s != "Phenotypic selection, 4-year cycle", ],
              aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se), alpha = 0.15, colour = NA) +
  geom_line(data = b[b$s != "Phenotypic selection, 4-year cycle", ], linewidth = 0.9) +
  geom_step(data = b[b$s == "Phenotypic selection, 4-year cycle", ], linewidth = 0.9) +
  annotate("segment", x = 4, xend = 4, y = 1.39, yend = 2.65, colour = "#C2410C", linewidth = 0.4,
           arrow = arrow(ends = "both", length = unit(2.5, "pt"))) +
  annotate("text", x = 4.4, y = 2.0, label = "year 4:\n2.6 vs 1.4", hjust = 0, size = 2.4, colour = "#C2410C",
           fontface = "bold", family = fam, lineheight = 0.9) +
  scale_colour_manual(values = c("#0F766E", "#334155")) + scale_fill_manual(values = c("#0F766E", "#334155"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Years", y = "Mean genetic value (base genetic SD)", title = "Early recycling raises gain sooner") + th

p <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = fam), plot.tag.position = c(0, 1.02))
ggsave("../figures/JXB_Fig4.pdf", p, width = 180, height = 82, units = "mm", device = cairo_pdf)
ggsave("../figures/JXB_Fig4.png", p, width = 180, height = 82, units = "mm", dpi = 600, type = "cairo")
