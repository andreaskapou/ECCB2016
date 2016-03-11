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
cgi_file <- "../datasets/ENCODE/hg19.cgi.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"

gm_rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-GM12878-rep1.bed"
h1_rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-H1-hESC-rep1.bed"
k562_rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-K562-rep1.bed"

chr_discarded <- c("chrX", "chrY", "chrM")
upstream      <- -1000
downstream    <- 1000

cgi_inters <- intersect_cgi_tss(cgi_file = cgi_file,
                                rna_files = k562_rnaseq_file,
                                chr_discarded = chr_discarded,
                                chrom_size_file = hg19_file,
                                upstream = upstream,
                                downstream = downstream)

# --------------------------------------
# Store the results
# --------------------------------------
filename <- paste0("../files/cgi_250_inter_ids", ".RData")
save(cgi_inters, file = filename)