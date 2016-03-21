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
source("init_regr_parameters.R")

rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsH1hescHaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-H1-hESC-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"


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
filename <- paste0("../files/model_corr_H1_5_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, out_prof_full, out_prof_rmse, out_prof_cpg, 
     out_prof, out_mean, file = filename)

# library(ggplot2)
# 
# # --------------------------------------
# # Compute Pearson's r
# # --------------------------------------
# prof_cor <- vector(mode = "numeric", length = iter)
# mean_cor <- vector(mode = "numeric", length = iter)
# for (i in 1:iter){
#   prof_cor[i] <- stats::cor(out_prof[[i]]$test_pred, out_prof[[i]]$test$y)
#   mean_cor[i] <- stats::cor(out_mean[[i]]$test_pred, out_mean[[i]]$test$y)
# }
# 
# corr_data <- data.frame(rbind(as.matrix(prof_cor), as.matrix(mean_cor)), 
#                         as.matrix(rep("H1-hESC", 2*iter)), 
#                         rbind(as.matrix(rep("Profile", iter)), as.matrix(rep("Mean", iter))),
#                         stringsAsFactors = TRUE)
# colnames(corr_data) <- c("R", "cell_line", "Method")
# 
# 
# ggplot(corr_data, aes(factor(Method), R)) + 
#   geom_boxplot() +
#   labs(list(x = "", y = "Correlation Coefficient")) + 
#   labs(title= "")

