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
k562_file <- "../files/cluster_K562_7000_5_4_FriMar041743.RData"
gm_file   <- "../files/cluster_GM_7000_5_4_FriMar041253.RData"
h1_file   <- "../files/cluster_H1_7000_5_4_FriMar041724.RData"

# -----------------------------------------
# Load saved data for K562
# -----------------------------------------
load(k562_file)
k562_HTS_data <- HTS_data
k562_proc_data <- proc_data
k562_mix_model <- mix_model

K <- k562_mix_model$K

# Cluster labels for K562 cell line
k562_labels <- list()
k562_expr <- list()
k562_gene_ids <- list()
for (i in 1:K){
  k562_labels[[i]] <- which(k562_mix_model$labels == i)
  k562_expr[[i]] <- k562_proc_data$Y[k562_labels[[i]]]
  k562_gene_ids[[i]] <- k562_proc_data$genes$ensembl_id[k562_labels[[i]]]
  #write(k562_gene_ids[[i]], paste0("../results/k562_7000_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}



# -----------------------------------------
# Load saved data for GM12878
# -----------------------------------------
load(gm_file)
gm_HTS_data <- HTS_data
gm_proc_data <- proc_data
gm_mix_model <- mix_model

# Cluster labels for K562 cell line
gm_labels <- list()
gm_expr <- list()
gm_gene_ids <- list()
for (i in 1:K){
  gm_labels[[i]] <- which(gm_mix_model$labels == i)
  gm_expr[[i]] <- gm_proc_data$Y[gm_labels[[i]]]
  gm_gene_ids[[i]] <- gm_proc_data$genes$ensembl_id[gm_labels[[i]]]
  #write(gm_gene_ids[[i]], paste0("../results/gm_7000_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}


# -----------------------------------------
# Load saved data for H1-hESC
# -----------------------------------------
load(h1_file)
h1_HTS_data <- HTS_data
h1_proc_data <- proc_data
h1_mix_model <- mix_model

# Cluster labels for K562 cell line
h1_labels <- list()
h1_expr <- list()
h1_gene_ids <- list()
for (i in 1:K){
  h1_labels[[i]] <- which(h1_mix_model$labels == i)
  h1_expr[[i]] <- h1_proc_data$Y[h1_labels[[i]]]
  h1_gene_ids[[i]] <- h1_proc_data$genes$ensembl_id[h1_labels[[i]]]
  #write(h1_gene_ids[[i]], paste0("../results/h1_7000_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}




# Plot methylation profiles for K = 5 clusters for DF Old mouse model
plot_cluster_prof(k562_mix_model, k562_mix_model$basis, FALSE)
# Corresponding gene expression levels for each cluster K
plot_cluster_box(k562_expr, FALSE)



# Plot methylation profiles for K = 5 clusters for DF Old mouse model
plot_cluster_prof(gm_mix_model, gm_mix_model$basis, FALSE)
# Corresponding gene expression levels for each cluster K
plot_cluster_box(gm_expr, FALSE)


# Plot methylation profiles for K = 5 clusters for DF Old mouse model
plot_cluster_prof(h1_mix_model, h1_mix_model$basis, FALSE)
# Corresponding gene expression levels for each cluster K
plot_cluster_box(h1_expr, FALSE)