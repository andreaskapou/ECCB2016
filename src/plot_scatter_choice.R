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
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)


# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
k562_file <- "../files/choice_corr_K562_ThuMar031032.RData"
load(k562_file)

k562_prof <- ggplot_scatt_regr_test2(out_prof, main_lab = "K562 Methylation Profile", is_margins = TRUE)
k562_mean <- ggplot_scatt_regr_test2(out_mean, main_lab = "K562 Mean Methylation", is_margins = TRUE)

corr_plot <- plot_grid(k562_prof, k562_mean, labels = c("A", "B"), 
                              label_size = 29, ncol = 2, nrow = 1)

save_plot("../figures/k562-scatter.pdf", corr_plot, ncol = 2, nrow = 1,
          base_height = 6.5, base_width = 6.5)

#ggsave("../figures/gg-k562-scatter-profile.pdf", k562_prof, width = 8, height = 8, units = "in")
#ggsave("../figures/gg-k562-scatter-mean.pdf", k562_mean, width = 8, height = 8, units = "in")

# t = 195
# plot(out_prof$basis, proc_data$obs[[t]], out_prof$W_opt[t,])
# plot(out_mean$basis, proc_data$obs[[t]], out_mean$W_opt[t,])

# plot_scatt_regr_test(out_prof, main_lab = "K562 Methylation Profile", is_margins = TRUE)
# plot_scatt_regr_test(out_mean, main_lab = "K562 Mean Methylation", is_margins = TRUE)

# plot_scatt_regr_train(out_prof, is_margins = TRUE)
# plot_scatt_regr_train(out_mean, main_lab = "Mean Methylation", is_margins = TRUE)