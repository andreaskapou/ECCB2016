# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(ggplot2)
library(cowplot)
library(gridExtra)
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)


# -----------------------------------------
# Initialize variables
# -----------------------------------------
k562_file <- "../files/model_corr_K562_ThuMar031216.RData"
gm_file   <- "../files/model_corr_GM_ThuMar031156.RData"
h1_file   <- "../files/model_corr_H1_ThuMar031152.RData"

# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
load(k562_file)

# Get the number of different dataset splits iterations
iter <- length(out_prof)

# Compute Pearson's r
k562_prof <- compute_corr(out       = out_prof, 
                          iter      = iter, 
                          cell_line = "K562",
                          method    = "Profile")
k562_prof_rmse <- compute_corr(out       = out_prof_rmse, 
                               iter      = iter, 
                               cell_line = "K562",
                               method    = "Profile+rmse")
k562_prof_cpg <- compute_corr(out       = out_prof_cpg, 
                              iter      = iter, 
                              cell_line = "K562",
                              method    = "Profile+#CpG")
k562_prof_full <- compute_corr(out       = out_prof_full, 
                               iter      = iter, 
                               cell_line = "K562",
                               method    = "Profile full")
k562_mean <- compute_corr(out       = out_mean, 
                          iter      = iter, 
                          cell_line = "K562",
                          method    = "Mean")

# -----------------------------------------
# Load saved data for GM12878
# -----------------------------------------
load(gm_file)

# Compute Pearson's r
gm_prof <- compute_corr(out       = out_prof, 
                        iter      = iter, 
                        cell_line = "GM12878",
                        method    = "Profile")
gm_prof_rmse <- compute_corr(out       = out_prof_rmse, 
                             iter      = iter, 
                             cell_line = "GM12878",
                             method    = "Profile+rmse")
gm_prof_cpg <- compute_corr(out       = out_prof_cpg, 
                            iter      = iter, 
                            cell_line = "GM12878",
                            method    = "Profile+#CpG")
gm_prof_full <- compute_corr(out       = out_prof_full, 
                             iter      = iter, 
                             cell_line = "GM12878",
                             method    = "Profile full")
gm_mean <- compute_corr(out       = out_mean, 
                        iter      = iter, 
                        cell_line = "GM12878",
                        method    = "Mean")


# -----------------------------------------
# Load saved data for H1-hESC
# -----------------------------------------
load(h1_file)

# Compute Pearson's r
h1_prof <- compute_corr(out       = out_prof, 
                        iter      = iter, 
                        cell_line = "H1-hESC",
                        method    = "Profile")
h1_prof_rmse <- compute_corr(out       = out_prof_rmse, 
                             iter      = iter, 
                             cell_line = "H1-hESC",
                             method    = "Profile+rmse")
h1_prof_cpg <- compute_corr(out       = out_prof_cpg, 
                            iter      = iter, 
                            cell_line = "H1-hESC",
                            method    = "Profile+#CpG")
h1_prof_full <- compute_corr(out       = out_prof_full, 
                             iter      = iter, 
                             cell_line = "H1-hESC",
                             method    = "Profile full")
h1_mean <- compute_corr(out       = out_mean, 
                        iter      = iter, 
                        cell_line = "H1-hESC",
                        method    = "Mean")


# ---------------------------------------
# Concatenate correlations from all cell lines
# and all methods
# ---------------------------------------
corr <- data.frame(rbind(k562_prof_full, k562_prof_rmse, k562_prof_cpg, 
                         k562_prof, k562_mean, 
                         gm_prof_full, gm_prof_rmse, gm_prof_cpg, 
                         gm_prof, gm_mean,
                         h1_prof_full, h1_prof_rmse, h1_prof_cpg, 
                         h1_prof, h1_mean))
colnames(corr) <- c("r", "cell_line", "method")
corr$method <- factor(corr$method, levels = c("Profile full", 
                                              "Profile+rmse", 
                                              "Profile+#CpG", 
                                              "Profile",
                                              "Mean"))


# ---------------------------------------
# Paired Wilcox test
# ---------------------------------------

# K562 against all other cell lines
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Profile full" & corr$cell_line == "GM12878"], 
            paired = TRUE, 
            alternative = "greater",
            exact = FALSE)$p.value  # 4.784587e-05
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Profile full" & corr$cell_line == "H1-hESC"], 
            paired = TRUE, 
            alternative = "greater",
            exact = FALSE)$p.value  # 4.784587e-05

# K562 profile against K562 mean
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Mean" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "greater",
            exact = FALSE)$p.value  # 4.784587e-05

# K562 profile full against K562 Profile
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Profile" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "greater",
            exact = FALSE)$p.value  # 4.784587e-05

# Profile against mean using all the cell lines
wilcox.test(corr$r[corr$method == "Profile full"], 
            corr$r[corr$method == "Mean"], 
            paired = TRUE, 
            alternative = "greater", 
            exact = FALSE)$p.value  # 8.356643e-12


# H1 against all other cell lines
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "H1-hESC"], 
            corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 4.784587e-05

wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "H1-hESC"], 
            corr$r[corr$method == "Profile full" & corr$cell_line == "GM12878"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 4.784587e-05




# K562 profile against K562 Profile+rmse
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Profile+rmse" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 4.784587e-05
# K562 profile against K562 Profile+#CpG
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Profile+#CpG" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 0.2220406



# GM12878 profile against GM12878 Profile+rmse
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "GM12878"], 
            corr$r[corr$method == "Profile+rmse" & corr$cell_line == "GM12878"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 4.784587e-05
# GM12878 profile against GM12878 Profile+#CpG
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "GM12878"], 
            corr$r[corr$method == "Profile+#CpG" & corr$cell_line == "GM12878"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 0.1801876



# H1-hESC profile against H1-hESC Profile+rmse
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "H1-hESC"], 
            corr$r[corr$method == "Profile+rmse" & corr$cell_line == "H1-hESC"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 0.0001363533
# H1-hESC profile against H1-hESC Profile+#CpG
wilcox.test(corr$r[corr$method == "Profile" & corr$cell_line == "H1-hESC"], 
            corr$r[corr$method == "Profile+#CpG" & corr$cell_line == "H1-hESC"], 
            paired = TRUE, 
            alternative = "less",
            exact = FALSE)$p.value  # 0.02094403


# ---------------------------------------
# Create box-plot comparing average methylation and methylation 
# profile, across all the  different cell lines.
# ---------------------------------------
model_boxplot <- ggplot(corr, aes(method, r)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(values=c("coral2","lightsalmon2","darkgoldenrod3",
                             "darkgoldenrod1","cornflowerblue"), 
                    name="Model",
                    #breaks=c("PrF", "PrR", "PrC", "Pr", "M"),
                    labels=c("Profile full", 
                             "Profile+rmse", 
                             "Profile+#CpG", 
                             "Profile",
                             "Mean")) +
  facet_grid( . ~ cell_line ) +
  theme_bw() +
  # coord_flip() + 
  labs(list(title= "", 
            x = "", 
            y = "Pearson's r")) +
  theme(axis.text.x = element_blank(), #element_text(size=17, angle=90, vjust = 0.4), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(face="bold", color = "black", size=23),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 20),
        text = element_text(size=21)) + 
  ggtitle("Model performances across cell lines")
  # guides(fill = FALSE)

# ggsave("../figures/model-corr-boxplot.pdf", model_boxplot, width = 11.5, height = 7, units = "in")




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# ------------------------------------
# Plot everything in the same figure
# ------------------------------------
k562_file <- "../files/choice_corr_K562_ThuMar031032.RData"

load(k562_file)
k562_prof <- ggplot_scatt_regr_test(out_prof, main_lab = "K562 Methylation Profile", is_margins = TRUE)
k562_mean <- ggplot_scatt_regr_test(out_mean, main_lab = "K562 Mean Methylation", is_margins = TRUE)

# grid.arrange(arrangeGrob(k562_prof, k562_mean), model_boxplot, ncol = 2)


corr_plot <- ggdraw() +
  draw_plot(k562_prof, 0, 0.5, 0.5, .5) +
  draw_plot(k562_mean, 0.5, 0.5, .5, .5) +
  draw_plot(model_boxplot, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "C", "B"), c(0, 0, 0.5), c(1, 0.5, 1), size = 22)

# save_plot("../figures/model-performance.pdf", corr_plot, ncol = 2, nrow = 2,
#           base_height = 6, base_width = 6)

