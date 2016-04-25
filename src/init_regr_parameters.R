# -------------------------------------------------
# Script for initializing experimental parameters
# for ENCODE datasets when performing regression.
# -------------------------------------------------


# -------------------------------------------------
# Initialize parameters for processing HTS data
# -------------------------------------------------
upstream    <- -7000
downstream  <- 7000
cpg_density <- 15
sd_thresh   <- 10e-02
min_bs_cov  <- 4
ignore_strand <- TRUE
chr_discarded <- c("chrX", "chrY", "chrM")

gene_expr_thresh <- FALSE
gene_outl_thresh <- TRUE
gene_log2_transf <- TRUE
is_fpkm <- TRUE
max_outl <- 600

# -------------------------------------------------
# Initialize parameters for regression model
# -------------------------------------------------
seed        <- 1234
formula     <- y ~ .
model_name  <- "svm"
train_perc  <- 0.7
opt_method  <- "CG"
opt_itnmax  <- 50
is_parallel <- TRUE
no_cores    <- 10
is_summary  <- TRUE

fit_feature   <- "RMSE"
cpg_dens_feat <- TRUE

basis_prof <- rbf.object(M = 9)
basis_mean <- polynomial.object(M = 0)