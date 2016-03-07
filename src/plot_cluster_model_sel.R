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
k562_file <- "../files/cluster_model_K562_5000_4_SatMar050033.RData"
gm_file   <- "../files/cluster_model_GM_5000_4_SatMar050011.RData"
h1_file   <- "../files/cluster_model_H1_5000_4_FriMar042317.RData"

# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
load(k562_file)
k562_HTS_data <- HTS_data
k562_proc_data <- proc_data
k562_mix_model <- mix_model

K <- k562_mix_model$K

k562_BIC <- vector(mode = "numeric", length = 9)
for (i in 2:10){
  k562_BIC[i-1] <- k562_mix_model[[i]]$BIC
}

plot(-k562_BIC,type="o", xlab="K", ylab="BIC score", col="red2", 
     main="K562 cell line", lwd=2)

# -----------------------------------------
# Load saved data for GM12878
# -----------------------------------------
load(gm_file)
gm_HTS_data <- HTS_data
gm_proc_data <- proc_data
gm_mix_model <- mix_model


gm_BIC <- vector(mode = "numeric", length = 9)
for (i in 2:10){
  gm_BIC[i-1] <- gm_mix_model[[i]]$BIC
}

plot(-gm_BIC,type="o", xlab="K", ylab="BIC score", col="red2", 
     main="GM12878 cell line", lwd=2)

# -----------------------------------------
# Load saved data for H1-hESC
# -----------------------------------------
load(h1_file)
h1_HTS_data <- HTS_data
h1_proc_data <- proc_data
h1_mix_model <- mix_model


h1_BIC <- vector(mode = "numeric", length = 9)
for (i in 2:10){
  h1_BIC[i-1] <- h1_mix_model[[i]]$BIC
}

plot(-h1_BIC,type="o", xlab="K", ylab="BIC score", col="red2", 
     main="H1-hESC cell line", lwd=2)
