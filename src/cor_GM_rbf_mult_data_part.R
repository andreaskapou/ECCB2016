# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)


# ------------------------------------------
# Initialize parameters
# ------------------------------------------
rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsGm12878HaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsGm12878HaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-GM12878-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"

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


# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
message("Promoter Region: ", downstream, "\n")
HTS_data <- process_haib_caltech(bs_files        = rrbs_file,
                                 rna_files       = rnaseq_file,
                                 chrom_size_file = hg19_file,
                                 chr_discarded   = chr_discarded,
                                 upstream        = upstream,
                                 downstream      = downstream,
                                 cpg_density     = cpg_density,
                                 sd_thresh       = sd_thresh,
                                 min_bs_cov      = min_bs_cov,
                                 ignore_strand   = ignore_strand)

proc_data <- preprocess_data(HTS_data         = HTS_data,
                             is_fpkm          = is_fpkm,
                             max_outl         = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)


# ------------------------------------------
# Final data:
#   X: contains methylation data for each promoter region
#   Y: contains the corresponding gene expression data
# ------------------------------------------
X <- proc_data$obs
Y <- proc_data$Y



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#--------------------------------------------
# Parameters for regression model
#--------------------------------------------
seed        <- 1234
formula     <- y ~ .
model_name  <- "svm"
train_perc  <- 0.7
opt_method  <- "CG"
opt_itnmax  <- 50
is_parallel <- TRUE
no_cores    <- 3
is_summary  <- TRUE

basis_prof <- rbf.object(M = 5, gamma = 14)
basis_mean <- polynomial.object(M = 0)

# Number of iterations
iter <- 20

# ------------------------------------------
# Apply methylation profiles full
# ------------------------------------------
set.seed(seed)
out_prof_full <- list()
for (i in 1:iter){
  out_prof_full[[i]] <- mpgex_regr(formula = formula,
                                   x = X,
                                   y = Y,
                                   model_name    = model_name,
                                   basis         = basis_prof,
                                   fit_feature   = "RMSE",
                                   cpg_dens_feat = TRUE,
                                   train_perc    = train_perc,
                                   opt_method    = opt_method,
                                   opt_itnmax    = opt_itnmax,
                                   is_parallel   = is_parallel,
                                   no_cores      = no_cores,
                                   is_summary    = is_summary)
}


# ------------------------------------------
# Apply methylation profiles RMSE
# ------------------------------------------
set.seed(seed)
out_prof_rmse <- list()
for (i in 1:iter){
  out_prof_rmse[[i]] <- mpgex_regr(formula = formula,
                                   x = X,
                                   y = Y,
                                   model_name    = model_name,
                                   basis         = basis_prof,
                                   fit_feature   = "RMSE",
                                   train_perc    = train_perc,
                                   opt_method    = opt_method,
                                   opt_itnmax    = opt_itnmax,
                                   is_parallel   = is_parallel,
                                   no_cores      = no_cores,
                                   is_summary    = is_summary)
}


# ------------------------------------------
# Apply methylation profiles CpG density
# ------------------------------------------
set.seed(seed)
out_prof_cpg <- list()
for (i in 1:iter){
  out_prof_cpg[[i]] <- mpgex_regr(formula = formula,
                                  x = X,
                                  y = Y,
                                  model_name    = model_name,
                                  basis         = basis_prof,
                                  cpg_dens_feat = TRUE,
                                  train_perc    = train_perc,
                                  opt_method    = opt_method,
                                  opt_itnmax    = opt_itnmax,
                                  is_parallel   = is_parallel,
                                  no_cores      = no_cores,
                                  is_summary    = is_summary)
}


# ------------------------------------------
# Apply methylation profiles 
# ------------------------------------------
set.seed(seed)
out_prof <- list()
for (i in 1:iter){
  out_prof[[i]] <- mpgex_regr(formula = formula,
                              x = X,
                              y = Y,
                              model_name    = model_name,
                              basis         = basis_prof,
                              train_perc    = train_perc,
                              opt_method    = opt_method,
                              opt_itnmax    = opt_itnmax,
                              is_parallel   = is_parallel,
                              no_cores      = no_cores,
                              is_summary    = is_summary)
}


# ------------------------------------------
# Apply mean methylation
# ------------------------------------------
set.seed(seed)
out_mean <- list()
for (i in 1:iter){
  out_mean[[i]] <- mpgex_regr(formula = formula,
                              x = X,
                              y = Y,
                              model_name    = model_name,
                              basis         = basis_mean,
                              train_perc    = train_perc,
                              opt_method    = opt_method,
                              opt_itnmax    = opt_itnmax,
                              is_parallel   = is_parallel,
                              no_cores      = no_cores,
                              is_summary    = is_summary)
}

# --------------------------------------
# Store the results
# --------------------------------------
filename <- paste0("../files/model_corr_GM_5_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, out_prof_full, out_prof_rmse, out_prof_cpg, 
     out_prof, out_mean, file = filename)

#
# # --------------------
# # Create plots
# t = 1569
#
# plot(out_prof$basis, X[[t]], out_prof$W_opt[t,])
# plot(out_mean$basis, X[[t]], out_mean$W_opt[t,])
#
# plot_scatt_regr_test(out_prof, main_lab = "GM12878 Methylation Profile", is_margins = TRUE)
# plot_scatt_regr_test(out_mean, main_lab = "GM12878 Mean Methylation", is_margins = TRUE)
#
# plot_scatt_regr_train(out_prof, is_margins = TRUE)
# plot_scatt_regr_train(out_mean, main_lab = "Mean Methylation", is_margins = TRUE)

