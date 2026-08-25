arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) stop("Expected results and output directories")

results_dir <- normalizePath(arguments[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(arguments[[2]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(scales))

ink <- "#161616"
text_gray <- "#404040"
grid_gray <- "#DEDAD3"
paper <- "#FFFFFF"
teal <- "#008B74"
orange <- "#D75B28"
purple <- "#7655B5"
gold <- "#D99A18"

prcc <- read.csv(file.path(results_dir, "model", "prcc.csv"), check.names = FALSE)
sobol <- read.csv(file.path(results_dir, "model", "sobol_indices.csv"), check.names = FALSE)
structural <- read.csv(file.path(results_dir, "model", "structural_contrasts.csv"), check.names = FALSE)

required_prcc <- c("parameter", "outcome", "prcc")
required_sobol <- c("parameter", "outcome", "ST")
required_structural <- c("model", "median_delta_fni")
if (!all(required_prcc %in% names(prcc))) stop("prcc.csv schema is invalid")
if (!all(required_sobol %in% names(sobol))) stop("sobol_indices.csv schema is invalid")
if (!all(required_structural %in% names(structural))) stop("structural_contrasts.csv schema is invalid")

fni <- prcc[prcc$outcome == "delta_fni", c("parameter", "prcc")]
names(fni)[2] <- "fni"
extent <- prcc[prcc$outcome == "delta_extent", c("parameter", "prcc")]
names(extent)[2] <- "extent"
paired <- merge(fni, extent, by = "parameter")
paired <- paired[order(abs(paired$fni), decreasing = TRUE), ][seq_len(min(12, nrow(paired))), ]
paired <- paired[order(paired$fni), ]
paired$parameter <- factor(paired$parameter, levels = paired$parameter)
long_prcc <- rbind(
  data.frame(parameter = paired$parameter, outcome = "Functional index", estimate = paired$fni),
  data.frame(parameter = paired$parameter, outcome = "Numerical extent", estimate = paired$extent)
)

base_theme <- theme_minimal(base_size = 10.5, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = paper, colour = NA),
    panel.background = element_rect(fill = paper, colour = NA),
    plot.title = element_text(colour = ink, face = "bold", size = 12.5, margin = margin(b = 8)),
    plot.tag = element_text(colour = ink, face = "bold", size = 12.5),
    axis.title = element_text(colour = ink, size = 9.5),
    axis.text = element_text(colour = text_gray, size = 8.4),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid_gray, linewidth = 0.35),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(colour = text_gray, size = 8.5),
    plot.margin = margin(8, 10, 8, 8)
  )

plot_prcc <- ggplot(long_prcc, aes(x = estimate, y = parameter, colour = outcome, shape = outcome)) +
  geom_vline(xintercept = 0, linewidth = 0.65, colour = ink) +
  geom_segment(aes(x = 0, xend = estimate, yend = parameter), linewidth = 0.85, alpha = 0.78) +
  geom_point(size = 2.7, stroke = 1.0) +
  scale_colour_manual(values = c("Functional index" = orange, "Numerical extent" = teal)) +
  scale_shape_manual(values = c("Functional index" = 16, "Numerical extent" = 0)) +
  scale_x_continuous(limits = c(-0.85, 0.85), breaks = seq(-0.8, 0.8, 0.4), labels = label_number(accuracy = 0.1, style_positive = "plus")) +
  labs(title = "A   Partial rank correlations", x = "PRCC", y = NULL) +
  base_theme

sobol_fni <- sobol[sobol$outcome == "delta_fni", ]
sobol_fni <- sobol_fni[order(sobol_fni$ST, decreasing = TRUE), ][seq_len(min(8, nrow(sobol_fni))), ]
sobol_fni <- sobol_fni[order(sobol_fni$ST), ]
sobol_fni$parameter <- factor(sobol_fni$parameter, levels = sobol_fni$parameter)

plot_sobol <- ggplot(sobol_fni, aes(x = ST, y = parameter)) +
  geom_segment(aes(x = 0, xend = ST, yend = parameter), linewidth = 4.8, colour = purple, alpha = 0.82, lineend = "round") +
  geom_point(size = 2.6, colour = ink) +
  geom_text(aes(label = sprintf("%.2f", ST)), hjust = 0, nudge_x = 0.025, colour = text_gray, size = 2.8) +
  scale_x_continuous(limits = c(0, max(sobol_fni$ST) * 1.32), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B   Total order sensitivity", x = "Sobol total order index", y = NULL) +
  base_theme +
  theme(legend.position = "none")

structure_order <- c("baseline", "proliferation", "maturation", "integration", "full")
structural <- structural[structural$model %in% structure_order, ]
structural$model <- factor(structural$model, levels = structure_order)
structural$fill <- ifelse(structural$model == "full", "Full model", "Reduced structure")

plot_structural <- ggplot(structural, aes(x = model, y = median_delta_fni, fill = fill)) +
  geom_col(width = 0.66, colour = paper, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", median_delta_fni)), vjust = -0.45, colour = ink, size = 3.0) +
  scale_fill_manual(values = c("Full model" = orange, "Reduced structure" = gold)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "C   Matched structural experiment", x = NULL, y = "Median delta FNI") +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = grid_gray, linewidth = 0.35),
    legend.position = "none"
  )

combined <- plot_prcc | (plot_sobol / plot_structural + plot_layout(heights = c(1.1, 0.9)))
combined <- combined +
  plot_layout(widths = c(1.35, 1.0)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = paper, colour = NA),
      plot.tag = element_text(colour = ink, face = "bold", size = 12.5)
    )
  )

png_path <- file.path(output_dir, "Fig3_sensitivity.png")
tif_path <- file.path(output_dir, "Fig3_sensitivity.tif")
ggsave(png_path, combined, width = 12, height = 7.5, units = "in", dpi = 600, bg = paper)
ggsave(tif_path, combined, width = 12, height = 7.5, units = "in", dpi = 600, bg = paper, compression = "lzw")
cat("Wrote ggplot2 sensitivity figure to", output_dir, "\n")
