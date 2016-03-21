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

rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsK562HaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-K562-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"

upstream    <- c(-1000, -2000, -3000, -4000, -5000, -6000, -7000, -8000, -9000, -10000)
downstream  <- c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

HTS_data <- list()
proc_data <- list()
out_prof <- list()
out_mean <- list()

for (i in 1:length(upstream)){
  # ------------------------------------------
  # Read and preprocess HTS files
  # ------------------------------------------
  message("Promoter Region: ", downstream, "\n")
  HTS_data[[i]] <- process_haib_caltech(bs_files        = rrbs_file,
                                        rna_files       = rnaseq_file,
                                        chrom_size_file = hg19_file,
                                        chr_discarded   = chr_discarded,
                                        upstream        = upstream[i],
                                        downstream      = downstream[i],
                                        cpg_density     = cpg_density,
                                        sd_thresh       = sd_thresh,
                                        min_bs_cov      = min_bs_cov,
                                        ignore_strand   = ignore_strand)
  
  proc_data[[i]] <- preprocess_data(HTS_data         = HTS_data[[i]],
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
  X <- proc_data[[i]]$obs
  Y <- proc_data[[i]]$Y



  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  

  # ------------------------------------------
  # Apply methylation profiles
  # ------------------------------------------
  set.seed(seed)
  out_prof[[i]] <- mpgex_regr(formula = formula,
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


  # ------------------------------------------
  # Apply mean methylation
  # ------------------------------------------
  set.seed(seed)
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
filename <- paste0("../files/region_K562_5_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, out_prof, out_mean, file = filename)
