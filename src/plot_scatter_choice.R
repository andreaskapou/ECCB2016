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

k562_prof <- ggplot_scatt_regr_test2(out_prof[[1]], main_lab = "K562 Methylation Profile", is_margins = TRUE)
k562_mean <- ggplot_scatt_regr_test2(out_mean, main_lab = "K562 Mean Methylation", is_margins = TRUE)

pvalue_prof = cor.test(as.vector(out_prof$test_pred), out_prof$test$y, alternative = "greater")$p.value
pvalue_mean = cor.test(as.vector(out_mean$test_pred), out_mean$test$y, alternative = "greater")$p.value

corr_plot <- plot_grid(k562_prof, k562_mean, labels = c("A", "B"), 
                              label_size = 29, ncol = 2, nrow = 1)

# save_plot("../figures/k562-scatter.pdf", corr_plot, ncol = 2, nrow = 1,
#           base_height = 6.5, base_width = 6.5)

# ggsave("../figures/gg-k562-scatter-profile.pdf", k562_prof, width = 8, height = 8, units = "in")
# ggsave("../figures/gg-k562-scatter-mean.pdf", k562_mean, width = 8, height = 8, units = "in")



# --------------------------------
# Statistical significance tests
# --------------------------------

k562_file <- "../files/corr_K562_WedMar021931.RData"
gm_file   <- "../files/corr_GM_WedMar021943.RData"
h1_file   <- "../files/corr_H1_WedMar021830.RData"

load(k562_file)
K562_profile <- out_prof[[1]]
K562_mean <- out_mean[[1]]

load(gm_file)
GM_profile <- out_prof[[1]]
GM_mean <- out_mean[[1]]

load(h1_file)
H1_profile <- out_prof[[1]]
H1_mean <- out_mean[[1]]


k562_prof = cor.test(K562_profile$test_pred, K562_profile$test$y, alternative = "greater")
k562_mean = cor.test(K562_mean$test_pred, K562_mean$test$y, alternative = "greater")

gm_prof = cor.test(GM_profile$test_pred, GM_profile$test$y, alternative = "greater")
gm_mean = cor.test(GM_mean$test_pred, GM_mean$test$y, alternative = "greater")

h1_prof = cor.test(H1_profile$test_pred, H1_profile$test$y, alternative = "greater")
h1_mean = cor.test(H1_mean$test_pred, H1_mean$test$y, alternative = "greater")