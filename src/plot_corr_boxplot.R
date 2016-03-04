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
k562_file <- "../files/corr_K562_WedMar021931.RData"
gm_file   <- "../files/corr_GM_WedMar021943.RData"
h1_file   <- "../files/corr_H1_WedMar021830.RData"


# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
load(k562_file)

# Get the number of different dataset splits iterations
iter <- length(out_prof)

# Compute Pearson's r
corr_k562 <- compute_corr_cell_line(out_prof, out_mean, iter, "K562")


# -----------------------------------------
# Load saved data for GM12878
# -----------------------------------------
load(gm_file)
# Compute Pearson's r
corr_gm <- compute_corr_cell_line(out_prof, out_mean, iter, "GM12878")


# -----------------------------------------
# Load saved data for H1-hESC
# -----------------------------------------
load(h1_file)
# Compute Pearson's r
corr_h1 <- compute_corr_cell_line(out_prof, out_mean, iter, "H1-hESC")


# ---------------------------------------
# Concatenate correlations from all cell lines
# ---------------------------------------
corr <- data.frame(rbind(corr_k562, corr_gm, corr_h1))
colnames(corr) <- c("r", "cell_line", "method")
corr$method <- factor(corr$method, levels = c("Methylation Profile", "Mean Methylation"))


# ---------------------------------------
# Create box-plot comparing average methylation and methylation 
# profile, across all the  different cell lines.
# ---------------------------------------
ggplot(corr, aes(method, r, fill = method)) +
  theme_bw() +
  geom_boxplot() + facet_grid(cell_line ~ . ) +
  labs(list(title= "", 
            x = "", 
            y = "Pearson's r between predicted and measured expression")) +
  theme(axis.text.x = element_text(size=17), 
        axis.text.y = element_text(size=12), text = element_text(size=17)) +
  guides(fill = FALSE)


# ---------------------------------------
# Paired Wilcox test
# ---------------------------------------
wilcox.test(corr$r[corr$method == "Methylation Profile" & corr$cell_line == "K562"], 
            corr$r[corr$method == "Mean Methylation" & corr$cell_line == "K562"], 
            paired = TRUE, 
            alternative = "greater", 
            exact = FALSE)$p.value  # 9.126857e-07

wilcox.test(corr$r[corr$method == "Methylation Profile"], 
            corr$r[corr$method == "Mean Methylation"], 
            paired = TRUE, 
            alternative = "greater", 
            exact = FALSE)$p.value  # 9.126857e-07