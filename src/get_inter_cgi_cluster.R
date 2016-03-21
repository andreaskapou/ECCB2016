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
# Read CGI file
# -----------------------------------------
cgi_file <- "../files/cgi_500_inter_ids.RData"
load(cgi_file)

# -----------------------------------------
# Read cluster files
# -----------------------------------------
k562_files <- c("../results/comm_k562_7000_5_4_clust_1.txt",
                "../results/comm_k562_7000_5_4_clust_2.txt",
                "../results/comm_k562_7000_5_4_clust_3.txt",
                "../results/comm_k562_7000_5_4_clust_4.txt",
                "../results/comm_k562_7000_5_4_clust_5.txt")

gm_files   <- c("../results/comm_gm_7000_5_4_clust_1.txt",
                "../results/comm_gm_7000_5_4_clust_2.txt",
                "../results/comm_gm_7000_5_4_clust_3.txt",
                "../results/comm_gm_7000_5_4_clust_4.txt",
                "../results/comm_gm_7000_5_4_clust_5.txt")

h1_files   <- c("../results/comm_h1_7000_5_4_clust_1.txt",
                "../results/comm_h1_7000_5_4_clust_2.txt",
                "../results/comm_h1_7000_5_4_clust_3.txt",
                "../results/comm_h1_7000_5_4_clust_4.txt",
                "../results/comm_h1_7000_5_4_clust_5.txt")

k562 <- list()
gm <- list()
h1 <- list()
for (i in 1:length(k562_files)){
  k562[[i]] <- read.table(file = k562_files[i])
  gm[[i]] <- read.table(file = gm_files[i])
  h1[[i]] <- read.table(file = h1_files[i])
}

inters_cell <- list()
for (i in 1:length(k562_files)){
  inters_cell[[i]] <- Reduce(intersect, list(k562[[i]]$V1, gm[[i]]$V1, h1[[i]]$V1))
  print(length(inters_cell[[i]]))
  message("Cluster ", i, " Venn intersection --> ", length(inters_cell[[i]])/length(gm[[i]]$V1))
}
message("\n\n")

cgi_clust_inters <- list()
for (i in 1:length(k562_files)){
  cgi_clust_inters[[i]] <- Reduce(intersect, list(cgi_inters$ensembl_id, inters_cell[[i]]))
  print(length(cgi_clust_inters[[i]]))
  message("Cluster ", i, " CGI intersection --> ", length(cgi_clust_inters[[i]])/length(inters_cell[[i]]))
}


