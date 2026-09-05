args <- commandArgs(trailingOnly = TRUE)
results <- args[1]
outdir <- args[2]
.libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
ink <- "#171717"
gray <- "#696969"
orange <- "#D9572B"
theme_study <- theme_minimal(base_size = 9) + theme(plot.title = element_text(size = 11, face = "bold", colour = ink), axis.title = element_text(colour = ink), axis.text = element_text(colour = gray), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#DDDDDD", linewidth = 0.3), strip.text = element_text(face = "bold", colour = ink), strip.background = element_blank(), legend.position = "bottom", legend.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))

shapley <- read.csv(file.path(results, "model", "channel_shapley_summary.csv"))
shapley$channel <- factor(shapley$channel, levels = rev(c("activation", "proliferation", "maturation", "integration", "efficacy", "survival")))
shapley$outcome <- factor(shapley$outcome, levels = c("delta_extent", "delta_fni"), labels = c("Numerical extent", "Functional index"))
p1 <- ggplot(shapley, aes(y = channel, x = mean, colour = outcome, shape = outcome)) + geom_vline(xintercept = 0, linewidth = 0.4, colour = ink) + geom_segment(aes(x = q025, xend = q975, yend = channel), position = position_dodge(width = 0.45), linewidth = 0.7) + geom_point(position = position_dodge(width = 0.45), size = 2.4, stroke = 0.8) + scale_colour_manual(values = c("Numerical extent" = gray, "Functional index" = orange)) + scale_shape_manual(values = c("Numerical extent" = 15, "Functional index" = 16)) + labs(title = "A  Exact channel attribution", x = "Shapley contribution", y = NULL) + theme_study

sobol <- read.csv(file.path(results, "model", "sobol_indices.csv"))
selected <- unique(unlist(lapply(split(sobol, sobol$outcome), function(x) head(x$parameter[order(x$ST, decreasing = TRUE)], 9))))
sobol <- sobol[sobol$parameter %in% selected, ]
sobol$parameter <- factor(sobol$parameter, levels = rev(unique(sobol$parameter[order(sobol$ST, decreasing = TRUE)])))
sobol$outcome <- factor(sobol$outcome, levels = c("delta_extent", "delta_fni"), labels = c("Numerical extent", "Functional index"))
parameter_names <- c(b_p = "Progenitor division", tx_prolif = "SSRI proliferation", d_p = "Progenitor loss", k_pn = "Neuroblast transition", tx_eff = "SSRI efficacy", eff_mean = "Mean cell efficacy", k_mg = "Cell integration", b_a = "Precursor division", tx_integration = "SSRI integration", k_qa = "Stem-cell activation", tx_maturation = "SSRI maturation", k_nm = "Neuron maturation", tx_survival = "SSRI survival")
p2 <- ggplot(sobol, aes(y = parameter, x = ST)) + geom_segment(aes(x = ST_low, xend = ST_high, yend = parameter), colour = gray, linewidth = 0.7) + geom_point(aes(colour = outcome), size = 2.3) + facet_wrap(~outcome, scales = "free_y") + scale_y_discrete(labels = parameter_names) + scale_colour_manual(values = c("Numerical extent" = gray, "Functional index" = orange), guide = "none") + labs(title = "B  Global total-order sensitivity", x = "Sobol total-order index", y = NULL) + theme_study + theme(axis.text.y = element_text(size = 7.5))

interactions <- read.csv(file.path(results, "model", "channel_pair_interaction_summary.csv"))
interactions <- interactions[interactions$outcome == "delta_fni", ]
levels_channels <- c("activation", "proliferation", "maturation", "integration", "efficacy", "survival")
interactions$channel_a <- factor(interactions$channel_a, levels = levels_channels)
interactions$channel_b <- factor(interactions$channel_b, levels = levels_channels)
p3 <- ggplot(interactions, aes(x = channel_b, y = channel_a, fill = mean)) + geom_tile(colour = "white", linewidth = 0.8) + geom_text(aes(label = sprintf("%.2f", mean)), size = 2.7, colour = ink) + scale_fill_gradient2(low = ink, mid = "white", high = orange, midpoint = 0, name = "Mean interaction") + labs(title = "C  Pairwise interactions for functional change", x = NULL, y = NULL) + theme_study + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

draw_all <- function() {
    grid::grid.newpage()
    layout <- grid::grid.layout(2, 2, heights = grid::unit(c(0.53, 0.47), "npc"))
    grid::pushViewport(grid::viewport(layout = layout))
    print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid::popViewport()
}

pdf(file.path(outdir, "Fig3_mechanism.pdf"), width = 12, height = 8.5, useDingbats = FALSE)
draw_all()
dev.off()
png(file.path(outdir, "Fig3_mechanism.png"), width = 7200, height = 5100, res = 600, bg = "white")
draw_all()
dev.off()
tiff(file.path(outdir, "Fig3_mechanism.tiff"), width = 7200, height = 5100, res = 600, compression = "lzw", bg = "white")
draw_all()
dev.off()
