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
                               method    = "Profile+RMSE")
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
                             method    = "Profile+RMSE")
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
                             method    = "Profile+RMSE")
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
                                              "Profile+RMSE", 
                                              "Profile+#CpG", 
                                              "Profile",
                                              "Mean"))


# ---------------------------------------
# Create box-plot comparing average methylation and methylation 
# profile, across all the  different cell lines.
# ---------------------------------------
model_boxplot <- ggplot(corr, aes(method, r)) +
  geom_boxplot(aes(fill = method)) + 
  #scale_fill_manual(values=c("blue","green","red","cyan",
  #                           "darkblue","darkgreen")) +
  facet_grid(cell_line ~ . ) +
  theme_bw() +
  labs(list(title= "", 
            x = "", 
            y = "Pearson's r between predicted and measured expression")) +
  theme(axis.text.x = element_text(size=10, angle=22, vjust = 0.7), 
        axis.text.y = element_text(size=10), 
        plot.title = element_text(face="bold", color = "black", size=15),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=14)) +
  guides(fill = FALSE) +
  ggtitle("Model performances across cell lines")

ggsave("../figures/model-corr-boxplot.pdf", model_boxplot, width = 7, height = 9.5, units = "in")

# ---------------------------------------
# Paired Wilcox test
# ---------------------------------------
wilcox.test(corr$r[corr$method == "Profile full" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Mean" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "greater", 
            exact = TRUE)$p.value  # 9.536743e-07

wilcox.test(corr$r[corr$method == "Profile full"], 
            corr$r[corr$method == "Mean"], 
            paired = TRUE, 
            alternative = "greater", 
            exact = TRUE)$p.value  # 8.673617e-19