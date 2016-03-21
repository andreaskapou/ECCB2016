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

#--------------------------------------------
# Parameters for regression model
#--------------------------------------------
seed        <- 1234
formula     <- y ~ .
model_name  <- "svm"
train_perc  <- 0.9
opt_method  <- "CG"
opt_itnmax  <- 50
is_parallel <- TRUE
no_cores    <- 10
is_summary  <- TRUE

basis_prof <- rbf.object(M = 5, gamma = 14)
basis_mean <- polynomial.object(M = 0)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# ------------------------------------------
# K562 cell line
# ------------------------------------------
rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsK562HaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-K562-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"


# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
message("Promoter Region: ", downstream, "\n")
K562_HTS_data <- process_haib_caltech(bs_files        = rrbs_file,
                                      rna_files       = rnaseq_file,
                                      chrom_size_file = hg19_file,
                                      chr_discarded   = chr_discarded,
                                      upstream        = upstream,
                                      downstream      = downstream,
                                      cpg_density     = cpg_density,
                                      sd_thresh       = sd_thresh,
                                      min_bs_cov      = min_bs_cov,
                                      ignore_strand   = ignore_strand)

K562_proc_data <- preprocess_data(HTS_data         = K562_HTS_data,
                                  is_fpkm          = is_fpkm,
                                  max_outl         = max_outl,
                                  gene_expr_thresh = gene_expr_thresh,
                                  gene_outl_thresh = gene_outl_thresh,
                                  gene_log2_transf = gene_log2_transf)

# ------------------------------------------
# Final preprocessed data
# ------------------------------------------
K562_X <- K562_proc_data$obs
K562_Y <- K562_proc_data$Y


# ------------------------------------------
# Apply methylation profiles
# ------------------------------------------
set.seed(seed)
K562_out_prof <- mpgex_regr(formula = formula,
                       x = K562_X,
                       y = K562_Y,
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
K562_out_mean <- mpgex_regr(formula = formula,
                       x = K562_X,
                       y = K562_Y,
                       model_name    = model_name,
                       basis         = basis_mean,
                       train_perc    = train_perc,
                       opt_method    = opt_method,
                       opt_itnmax    = opt_itnmax,
                       is_parallel   = is_parallel,
                       no_cores      = no_cores,
                       is_summary    = is_summary)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# ------------------------------------------
# GM12878 cell line
# ------------------------------------------
rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsGm12878HaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsGm12878HaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-GM12878-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"


# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
message("Promoter Region: ", downstream, "\n")
GM_HTS_data <- process_haib_caltech(bs_files        = rrbs_file,
                                    rna_files       = rnaseq_file,
                                    chrom_size_file = hg19_file,
                                    chr_discarded   = chr_discarded,
                                    upstream        = upstream,
                                    downstream      = downstream,
                                    cpg_density     = cpg_density,
                                    sd_thresh       = sd_thresh,
                                    min_bs_cov      = min_bs_cov,
                                    ignore_strand   = ignore_strand)

GM_proc_data <- preprocess_data(HTS_data         = GM_HTS_data,
                                is_fpkm          = is_fpkm,
                                max_outl         = max_outl,
                                gene_expr_thresh = gene_expr_thresh,
                                gene_outl_thresh = gene_outl_thresh,
                                gene_log2_transf = gene_log2_transf)

# ------------------------------------------
# Final preprocessed data
# ------------------------------------------
GM_X <- GM_proc_data$obs
GM_Y <- GM_proc_data$Y


# ------------------------------------------
# Apply methylation profiles
# ------------------------------------------
set.seed(seed)
GM_out_prof <- mpgex_regr(formula = formula,
                          x = GM_X,
                          y = GM_Y,
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
GM_out_mean <- mpgex_regr(formula = formula,
                          x = GM_X,
                          y = GM_Y,
                          model_name    = model_name,
                          basis         = basis_mean,
                          train_perc    = train_perc,
                          opt_method    = opt_method,
                          opt_itnmax    = opt_itnmax,
                          is_parallel   = is_parallel,
                          no_cores      = no_cores,
                          is_summary    = is_summary)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# ------------------------------------------
# H1-hESC cell line
# ------------------------------------------
rrbs_file   <- c("../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz",
                 "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsH1hescHaibSitesRep2.bed.gz")
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-H1-hESC-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"


# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
message("Promoter Region: ", downstream, "\n")
H1_HTS_data <- process_haib_caltech(bs_files        = rrbs_file,
                                    rna_files       = rnaseq_file,
                                    chrom_size_file = hg19_file,
                                    chr_discarded   = chr_discarded,
                                    upstream        = upstream,
                                    downstream      = downstream,
                                    cpg_density     = cpg_density,
                                    sd_thresh       = sd_thresh,
                                    min_bs_cov      = min_bs_cov,
                                    ignore_strand   = ignore_strand)

H1_proc_data <- preprocess_data(HTS_data         = H1_HTS_data,
                                is_fpkm          = is_fpkm,
                                max_outl         = max_outl,
                                gene_expr_thresh = gene_expr_thresh,
                                gene_outl_thresh = gene_outl_thresh,
                                gene_log2_transf = gene_log2_transf)

# ------------------------------------------
# Final preprocessed data
# ------------------------------------------
H1_X <- H1_proc_data$obs
H1_Y <- H1_proc_data$Y


# ------------------------------------------
# Apply methylation profiles
# ------------------------------------------
set.seed(seed)
H1_out_prof <- mpgex_regr(formula = formula,
                          x = H1_X,
                          y = H1_Y,
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
H1_out_mean <- mpgex_regr(formula = formula,
                          x = H1_X,
                          y = H1_Y,
                          model_name    = model_name,
                          basis         = basis_mean,
                          train_perc    = train_perc,
                          opt_method    = opt_method,
                          opt_itnmax    = opt_itnmax,
                          is_parallel   = is_parallel,
                          no_cores      = no_cores,
                          is_summary    = is_summary)


# --------------------------------------
# Store the results
# --------------------------------------
filename <- paste0("../files/across_cell_line_corr_5_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(K562_HTS_data, K562_proc_data, K562_out_prof, K562_out_mean, 
     GM_HTS_data, GM_proc_data, GM_out_prof, GM_out_mean, 
     H1_HTS_data, H1_proc_data, H1_out_prof, H1_out_mean,
     file = filename)
