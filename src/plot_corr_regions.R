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
k562_file <- "../files/region_K562_WedMar022118.RData"
gm_file   <- "../files/region_GM_WedMar022115.RData"
h1_file   <- "../files/region_H1_WedMar022115.RData"


# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
load(k562_file)
# Get the number of different regions
iter <- length(out_prof)

# Compute Pearson's r
k562_prof <- compute_corr(out       = out_prof, 
                          iter      = iter, 
                          cell_line = "K562",
                          method    = "Profile")

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
h1_mean <- compute_corr(out       = out_mean, 
                        iter      = iter, 
                        cell_line = "H1-hESC",
                        method    = "Mean")