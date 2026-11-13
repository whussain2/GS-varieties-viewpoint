# Counting varieties misjudges genomic selection: what it delivers, and when

R and Python source code and simulated data for the Viewpoint by **Waseem Hussain** and **Vinay Kumar Nandigam** (PlantNura Consulting & Solutions LLP), submitted to *Journal of Experimental Botany*.

> Genomic selection (GS) does not make varieties, and it is not meant to. Its job is to discard poor lines early, using prediction from a mixed model that borrows information from relatives through marker-based covariance, and to recycle better parents faster. Its product is a better breeding population, which supplies product development, where varieties are identified in the field.

## Contents

| Script | Reproduces | Output |
|---|---|---|
| `code/Fig1_Fig2_simulate_program.R` | One simulated program, start to finish (Box 1; Figs 1 and 2): crosses → F5 lines → SNPs → **G** matrix → CDmean training set → multi-location trial → RR-BLUP, GBLUP, BayesB → ranked GEBVs | `Fig1_data.json`, `Fig1_*.csv` |
| `code/Box1_100_replicate_programmes.R` | The 100-program statistics quoted in the text (shrinkage, reliability, realized accuracy across and within families, selection outcomes) | `Box1_100_replicates*.csv` |
| `code/Fig3C_allele_frequency_sim.R` | Favorable-allele frequencies under recurrent selection (Fig. 3C) | `Fig3C_*.csv` |
| `code/Fig4A_selection_curves.R` | Selection on a criterion of accuracy *r*: a good set is kept, the best line rarely (Fig. 4A) | `Fig4A_selection_curves.csv` |
| `code/Fig4B_multicycle_GS_vs_PS.R` | 20-year simulation, GS with early recycling vs phenotypic selection, 30 replicates (Fig. 4B) | `Fig4B_multicycle_summary.csv` |
| `code/Fig4_plot.R` | Draws Fig. 4 | `figures/JXB_Fig4.pdf/.png` |
| `figures_src/build_fig1.py`, `build_fig2.py`, `render.py` | Compose Figs 1–3 as HTML/SVG from the simulated data and render them to PDF/PNG | `figures/JXB_Fig1–3` |

All simulations are seeded and written in plain base R so that each step can be read and followed, with **rrBLUP** for REML and RR-BLUP, **BGLR** for BayesB, and **ggplot2**/**patchwork** for Fig. 4. Only the summary tables behind the numbers quoted in the paper are stored in `code/`; every other output (genotypes, **G** matrix, field trial, GEBVs, per-replicate results, allele frequencies) is regenerated exactly by the seeded scripts. Run the R scripts before the Python figure scripts.

## How to reproduce

```r
# R 4.3.3; install.packages(c("rrBLUP", "BGLR", "jsonlite", "ggplot2", "patchwork"))
setwd("code")
source("Fig1_Fig2_simulate_program.R")
source("Box1_100_replicate_programmes.R")
source("Fig3C_allele_frequency_sim.R")
source("Fig4A_selection_curves.R")
source("Fig4B_multicycle_GS_vs_PS.R")
source("Fig4_plot.R")
```

Figures 1–3 are composed as HTML from the simulated data and rendered with Python and Playwright (Chromium):

```bash
cd figures_src
python3 build_fig1.py && python3 build_fig2.py
python3 render.py "$PWD/Fig1_composite.html" "$PWD/../figures/JXB_Fig1"
python3 render.py "$PWD/Fig2_composite.html" "$PWD/../figures/JXB_Fig2"
python3 render.py "$PWD/Fig3_composite.html" "$PWD/../figures/JXB_Fig3"
```

Fig. 3B probability: the gamete from Plant 2 = (1 − 0.02)(1 − 0.40)/2 × 0.40(1 − 0.40)/2 ≈ 0.035; the gamete from Plant 4 = (1 − 0.02)(0.40)/2 × 0.40(1 − 0.40)/2 ≈ 0.024; product ≈ 0.00083, about 1 in 1,200.

## Software versions

R 4.3.3 · rrBLUP 4.6.3 · BGLR 1.1.4 · ggplot2 4.0.3 · Python 3 with Playwright

## Citation

Hussain W, Nandigam VK. Counting varieties misjudges genomic selection: what it delivers, and when. *Journal of Experimental Botany* (submitted). Code: https://github.com/whussain2/GS-varieties-viewpoint

## License

Code is released under the MIT License (see `LICENSE`). Simulated data may be reused with attribution.
