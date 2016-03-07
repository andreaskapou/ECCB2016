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

upstream    <- -5000
downstream  <- 5000
cpg_density <- 15
sd_thresh   <- 7e-02
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


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#--------------------------------------------
# Parameters for EM algorithm
#--------------------------------------------
seed        <- 1234
K           <- 5
pi_k        <- NULL
w           <- NULL
basis       <- polynomial.object(M = 4)
em_max_iter <- 20
epsilon_conv <- 1e-4
opt_method  <- "CG"
opt_itnmax  <- 50
init_opt_itnmax <- 50
is_parallel <- TRUE
no_cores    <- 5
is_verbose  <- TRUE

mix_model <- mpgex_cluster(x     = proc_data$obs,
                           K     = K,
                           pi_k  = pi_k,
                           w     = w,
                           basis = basis,
                           em_max_iter  = em_max_iter,
                           epsilon_conv = epsilon_conv,
                           opt_method   = opt_method,
                           opt_itnmax   = opt_itnmax,
                           init_opt_itnmax = init_opt_itnmax,
                           is_parallel  = is_parallel,
                           no_cores     = no_cores,
                           is_verbose   = is_verbose)


# --------------------------------------
# Store the results
# --------------------------------------
filename <- paste0("../files/cluster_poly_GM_5000_5_5_7_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, mix_model, file = filename)
