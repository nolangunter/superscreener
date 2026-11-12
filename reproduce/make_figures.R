# make_figures.R
# ---------------------------------------------------------------------------
# Simulation figures of the chapter, from the CSVs written by summarize_simulation.R.
#   fig1_sim_prediction_by_K  AUC shortfall at a fixed K: SuperScreener's choice of ranking vs consensus
#   fig2_sim_stability_by_K   stability at a fixed K: SuperScreener, consensus, univariate LRT
#   fig3_sim_choosing_K       variables used and AUC shortfall against the price c
# Usage (from the repository root): Rscript reproduce/make_figures.R [--outdir reproduce/output/sim]
# Needs dplyr, tidyr, readr and ggplot2.
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
SUMDIR <- file.path(path.expand(arg("--outdir", "reproduce/output/sim")), "summary")
OUT <- file.path(SUMDIR, "figures"); dir.create(OUT, showWarnings = FALSE)

ORANGE <- "#eb6834"; BLUE <- "#2a78d6"; AQUA <- "#1baf7a"; INK <- "#30343b"; MUTED <- "#5b6370"; GRAY <- "#8c929b"
theme_c <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#e6e8eb", linewidth = 0.3),
        axis.text = element_text(color = MUTED), axis.title = element_text(color = INK),
        strip.text = element_text(face = "bold", hjust = 0, size = 10.5), legend.position = "top",
        legend.title = element_blank(), panel.spacing = unit(1.2, "lines"))
save <- function(g, name, w = 8, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), g, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), g, width = w, height = h, dpi = 220, bg = "white")
}
SCEN <- c(n500 = "Moderate sample (n = 500, p = 300)", prox_high = "High proxy correlation (n = 1000, p = 300)",
          small_n = "Very small sample (n = 100, p = 2000)", strong_nl = "Strongly nonlinear (n = 200, p = 500)")
scen_f <- function(x) factor(SCEN[x], SCEN[names(SCEN) %in% x])
xlab_K <- "Number of columns kept (K), fixed by the analyst"

## Figure 1
pred <- read_csv(file.path(SUMDIR, "pred_all.csv"), show_col_types = FALSE)
K_GRID <- sort(unique(pred$K[pred$method == "consensus"]))
p1 <- pred |> filter(method %in% c("consensus", "ss"), !is.na(K), K > 0) |>
  mutate(Method = factor(ifelse(method == "ss", "SuperScreener", "Consensus ranking"), c("SuperScreener", "Consensus ranking")),
         scenario = scen_f(scenario))
full1 <- pred |> filter(method == "full") |> mutate(scenario = scen_f(scenario))
g1 <- ggplot(p1, aes(K, short_m, color = Method)) +
  geom_hline(data = full1, aes(yintercept = short_m, linetype = "Full model"), color = GRAY, linewidth = 0.6) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  facet_wrap(~scenario, nrow = 2, scales = "free_y") +
  scale_color_manual(values = c(ORANGE, BLUE)) + scale_linetype_manual(values = "dashed") +
  scale_x_continuous(breaks = K_GRID) +
  labs(x = xlab_K, y = "AUC shortfall from the true risk (x 100)") + theme_c
save(g1, "fig1_sim_prediction_by_K")

## Figure 2
stab <- read_csv(file.path(SUMDIR, "stab_all.csv"), show_col_types = FALSE)
LV2 <- c(ss = "SuperScreener", consensus = "Consensus ranking", lrt = "Univariate LRT ranking")
p2 <- stab |> filter(method %in% names(LV2), K > 0) |>
  mutate(Method = factor(LV2[method], LV2), scenario = scen_f(scenario))
g2 <- ggplot(p2, aes(K, phi_m, color = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) + facet_wrap(~scenario, nrow = 2) +
  scale_color_manual(values = c(ORANGE, BLUE, AQUA)) + scale_x_continuous(breaks = K_GRID) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = xlab_K, y = "Stability of the selected set (Nogueira index)") + theme_c
save(g2, "fig2_sim_stability_by_K")

## Figure 3: the price on the x axis (c = 0 has no place on a log axis); the applied rule and the
## full model do not depend on the price, so they are flat reference lines
ks <- read_csv(file.path(SUMDIR, "kselect_summary.csv"), show_col_types = FALSE) |> mutate(scenario = scen_f(scenario))
LV3 <- c("SuperScreener", "SuperScreener, at least 5 variables", "Applied rule", "Full model")
MEAS <- c("Variables used", "AUC shortfall (x 100)")
long3 <- function(d) d |> transmute(scenario, c, Line, v1 = vars_m, v2 = short_m) |>
  pivot_longer(c(v1, v2), names_to = "measure") |> mutate(measure = factor(MEAS[match(measure, c("v1", "v2"))], MEAS))
path <- ks |> filter(rule %in% c("ss_full", "ss_floor"), c > 0) |>
  mutate(Line = factor(ifelse(rule == "ss_full", LV3[1], LV3[2]), LV3)) |> long3()
refs <- ks |> filter(rule %in% c("applied", "full"), abs(c - 0.01) < 1e-9) |>
  mutate(Line = factor(ifelse(rule == "applied", LV3[3], LV3[4]), LV3)) |> long3() |>
  filter(!(Line == "Full model" & measure == "Variables used"))
COL3 <- setNames(c(ORANGE, "#f4ab8c", BLUE, GRAY), LV3)
g3 <- ggplot() +
  geom_vline(xintercept = 1, color = "#d5d8dc", linewidth = 0.4) +
  geom_hline(data = refs, aes(yintercept = value, color = Line, linetype = Line), linewidth = 0.7) +
  geom_line(data = path, aes(100 * c, value, color = Line, linetype = Line), linewidth = 1) +
  geom_point(data = path, aes(100 * c, value, color = Line), size = 1.8, show.legend = FALSE) +
  facet_grid(measure ~ scenario, scales = "free_y", switch = "y",
             labeller = labeller(scenario = function(x) sub(" \\(", "\n(", x))) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  scale_color_manual(values = COL3, breaks = LV3) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "dashed"), LV3), breaks = LV3) +
  labs(x = "Price c (% of the best candidate's risk reduction, log scale)", y = NULL) + theme_c +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5, face = "plain"),
        strip.text.x = element_text(size = 9.5))
save(g3, "fig3_sim_choosing_K", w = 10, h = 5.6)
message("figures in ", OUT)
