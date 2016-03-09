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
k562_file <- "../files/cluster_K562_7000_5_3_7_MonMar071952.RData"
gm_file   <- "../files/cluster_GM_7000_5_3_7_MonMar071942.RData"
h1_file   <- "../files/cluster_H1_7000_5_3_7_MonMar071944.RData"

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
  print(length(k562_labels[[i]]))
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
  print(length(gm_labels[[i]]))
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
  print(length(h1_labels[[i]]))
  h1_expr[[i]] <- h1_proc_data$Y[h1_labels[[i]]]
  h1_gene_ids[[i]] <- h1_proc_data$genes$ensembl_id[h1_labels[[i]]]
  #write(h1_gene_ids[[i]], paste0("../results/h1_7000_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -------------------------------------
# Preprocess data for plotting
# -------------------------------------
merged_meth <- list(k562_mix_model, gm_mix_model, h1_mix_model)
merged_expr <- list(k562_expr, gm_expr, h1_expr)
cell_lines <- c("K562", "GM12878", "H1-hESC")
x_len   <- 2000


# ----------------------------------------------
# Create data frame containing all experimental 
# output for clustered methylation profiles
# ----------------------------------------------
df_meth <- data.frame(xs = numeric(),
                      y = numeric(), 
                      cluster = character(), 
                      cell_line = character(),
                      stringsAsFactors=TRUE)
for (i in 1:length(merged_meth)){
  for (k in 1:K){
    xs <- seq(-1, 1,len = x_len) # create some values
    y <- as.vector(eval_probit_function(merged_meth[[i]]$basis, 
                                        xs, 
                                        merged_meth[[i]]$w[, k]))
    cluster <- paste("Cluster", k)
    cell_line <- cell_lines[i]
    dd <- data.frame(xs, y, cluster, cell_line, stringsAsFactors = TRUE)
    df_meth <- data.frame(rbind(df_meth, dd))
  }
}


# ----------------------------------------------
# Create data frame containing all experimental 
# output for  clustered gene expression levels
# ----------------------------------------------
df_expr <- data.frame(expr = numeric(), 
                      cluster = character(), 
                      cell_line = character(),
                      stringsAsFactors=TRUE)

for (i in 1:length(merged_expr)){
  for (k in 1:K){
    expr <- merged_expr[[i]][[k]]
    cluster <- paste("Cluster", k)
    cell_line <- cell_lines[i]
    dd <- data.frame(expr, cluster, cell_line, stringsAsFactors = TRUE)
    df_expr <- data.frame(rbind(df_expr, dd))
  }
}



# ---------------------------------------------
# Create plots
# ---------------------------------------------
gg_prof <- ggplot_cluster_prof(df = df_meth)

gg_expr <- ggplot_cluster_expr(df = df_expr)