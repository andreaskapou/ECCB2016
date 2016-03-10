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
#----------------
# Create Venn diagrams for K = 5, RBF = 4
#----------------

# Cluster 1 --------------
k562_1 <- "../results/comm_k562_7000_5_4_clust_1.txt"
gm_1   <- "../results/comm_gm_7000_5_4_clust_1.txt"
h1_1   <- "../results/comm_h1_7000_5_4_clust_1.txt"

create_triple_venn(k562_1, gm_1, h1_1,
                   filename = "../figures/comm_cluster_1_7000_5_4.png")


# Cluster 2 --------------------
k562_2 <- "../results/comm_k562_7000_5_4_clust_2.txt"
gm_2   <- "../results/comm_gm_7000_5_4_clust_2.txt"
h1_2   <- "../results/comm_h1_7000_5_4_clust_2.txt"

create_triple_venn(k562_2, gm_2, h1_2,
                   filename = "../figures/comm_cluster_2_7000_5_4.png")



# Cluster 3 -------------------
k562_3 <- "../results/comm_k562_7000_5_4_clust_3.txt"
gm_3   <- "../results/comm_gm_7000_5_4_clust_3.txt"
h1_3   <- "../results/comm_h1_7000_5_4_clust_3.txt"

create_triple_venn(k562_3, gm_3, h1_3,
                   filename = "../figures/comm_cluster_3_7000_5_4.png")


# Cluster 4 -------------------
k562_4 <- "../results/comm_k562_7000_5_4_clust_4.txt"
gm_4   <- "../results/comm_gm_7000_5_4_clust_4.txt"
h1_4   <- "../results/comm_h1_7000_5_4_clust_4.txt"

create_triple_venn(k562_4, gm_4, h1_4,
                   filename = "../figures/comm_cluster_4_7000_5_4.png")


# Cluster 5 -------------------
k562_5 <- "../results/comm_k562_7000_5_4_clust_5.txt"
gm_5   <- "../results/comm_gm_7000_5_4_clust_5.txt"
h1_5   <- "../results/comm_h1_7000_5_4_clust_5.txt"

create_triple_venn(k562_5, gm_5, h1_5,
                   filename = "../figures/comm_cluster_5_7000_5_4.png")
