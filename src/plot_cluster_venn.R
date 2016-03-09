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
k562_1 <- "../results/k562_7000_5_4_clust_1_WedMar091917.txt"
gm_1   <- "../results/gm_7000_5_4_clust_1_WedMar091917.txt"
h1_1   <- "../results/h1_7000_5_4_clust_1_WedMar091917.txt"

create_triple_venn(k562_1, gm_1, h1_1,
                 filename = "../figures/Cluster_1_7000_5_4.png")


# Cluster 2 --------------------
k562_2 <- "../results/k562_7000_5_4_clust_2_WedMar091917.txt"
gm_2   <- "../results/gm_7000_5_4_clust_2_WedMar091917.txt"
h1_2   <- "../results/h1_7000_5_4_clust_2_WedMar091917.txt"

create_triple_venn(k562_2, gm_2, h1_2,
                   filename = "../figures/Cluster_2_7000_5_4.png")



# Cluster 3 -------------------
k562_3 <- "../results/k562_7000_5_4_clust_3_WedMar091917.txt"
gm_3   <- "../results/gm_7000_5_4_clust_3_WedMar091917.txt"
h1_3   <- "../results/h1_7000_5_4_clust_3_WedMar091917.txt"

create_triple_venn(k562_3, gm_3, h1_3,
                   filename = "../figures/Cluster_3_7000_5_4.png")


# Cluster 4 -------------------
k562_4 <- "../results/k562_7000_5_4_clust_4_WedMar091917.txt"
gm_4   <- "../results/gm_7000_5_4_clust_4_WedMar091917.txt"
h1_4   <- "../results/h1_7000_5_4_clust_4_WedMar091917.txt"

create_triple_venn(k562_4, gm_4, h1_4,
                   filename = "../figures/Cluster_4_7000_5_4.png")


# Cluster 5 -------------------
k562_5 <- "../results/k562_7000_5_4_clust_5_WedMar091917.txt"
gm_5   <- "../results/gm_7000_5_4_clust_5_WedMar091917.txt"
h1_5   <- "../results/h1_7000_5_4_clust_5_WedMar091917.txt"

create_triple_venn(k562_5, gm_5, h1_5,
                   filename = "../figures/Cluster_5_7000_5_4.png")