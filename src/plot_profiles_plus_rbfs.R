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
rrbs_file   <- "../datasets/ENCODE/BS-Seq/wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed.gz"
rnaseq_file <- "../datasets/ENCODE/RNA-Seq/GENCODE-v3-K562-rep1.bed"
hg19_file   <- "../datasets/ENCODE/hg19.chrom.sizes"

upstream    <- -7000
downstream  <- 7000
cpg_density <- 40
sd_thresh   <- 3e-02
min_bs_cov  <- 1
ignore_strand <- FALSE
chr_discarded <- c("chrX", "chrY", "chrM")

gene_expr_thresh <- FALSE
gene_outl_thresh <- FALSE
gene_log2_transf <- TRUE
is_fpkm <- FALSE
max_outl <- 600


# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
message("Promoter Region: ", downstream, "\n")
HTS_data <- process_haib_caltech_wrap(bs_files        = rrbs_file,
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
seed        <- 12345
formula     <- y ~ .
model_name  <- "svm"
train_perc  <- 0.6
opt_method  <- "CG"
opt_itnmax  <- 100
is_parallel <- TRUE
no_cores    <- 8
is_summary  <- TRUE

basis_prof <- rbf.object(M = 4, gamma = 12)
basis_mean <- polynomial.object(M = 0)


# ------------------------------------------
# Apply methylation profiles full
# ------------------------------------------
set.seed(seed)
out_prof <- mpgex_regr(formula = formula,
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
out_mean <- mpgex_regr(formula = formula,
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



# ------------------------------------------
# Create plots for FAM20C and DOCK5 gene promoters
# ------------------------------------------
xs <- seq(-1,1,len=2000) # Create some values

#
# CCR10 methylation profile
#
par(cex=1.25, mai=c(1.37,1.37,.7,.3) )
t <- 1465
proc_data$genes$gene_name[t]
x <- HTS_data$methyl_region[[t]][,1]
y <- HTS_data$methyl_region[[t]][,3]/HTS_data$methyl_region[[t]][,2]
xl <- seq(-1, 1, 0.001)
plot(x, y, col = "blue3", pch = 22, ylim = c(0,1), xlim = c(-1,1), lwd = 2,
     xlab = NA, ylab = NA, cex.axis = 1.3, xaxt = "n")
mtext(side = 1, "genomic region", line = 3, cex = 1.8)
mtext(side = 2, "methylation level", line = 3, cex = 1.8)
axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
title(main="Gene CCR10", line = 1, cex.main=1.6)
lines(x = xl, y = eval_probit_function(out_mean$basis, xl, out_mean$W_opt[t,]), 
      col = 'cornflowerblue', lwd = 3, lty = 2)
lines(x = xl, y = eval_probit_function(out_prof$basis, xl, out_prof$W_opt[t,]), 
      col = 'blue2', lwd = 3)

# dev.print(pdf, "../figures/ccr10-2.pdf", width = 10.5, height = 7)
# dev.off()



plot(NULL, xlim = c(-1,1), ylim = c(0,1), type = "n")
points(basis_prof$mus, rep(0, basis_prof$M), col=1:basis_prof$M, pch=19)
x.coord = seq(-1, 1, length = 250)
for (i in 1:basis_prof$M) {
  points(x.coord, qnorm(pnorm(apply(X = as.matrix(x.coord), 
                        MARGIN = 1,
                        FUN = rbf_basis,
                        mus = basis_prof$mus[i],
                        gamma = basis_prof$gamma))), 
         type="l", col=i)
}


plot(NULL, xlim = c(-1,1), ylim = c(0,1), type = "n")
points(basis_prof$mus, rep(0, basis_prof$M), col=1:basis_prof$M, pch=19)
x.coord = seq(-1, 1, length = 250)
for (i in 1:basis_prof$M) {
  evals <- pnorm(out_prof$W_opt[t,1] * apply(X = as.matrix(x.coord), 
                                               MARGIN = 1,
                                               FUN = rbf_basis,
                                               mus = basis_prof$mus[i],
                                               gamma = basis_prof$gamma) + 
                   out_prof$W_opt[t,i+1] * apply(X = as.matrix(x.coord), 
                                               MARGIN = 1,
                                               FUN = rbf_basis,
                                               mus = basis_prof$mus[i],
                                               gamma = basis_prof$gamma))
  points(x.coord, evals, 
         type="l", col=i)
}



#
# PLEKHH3 methylation profile
#
par(mfrow = c(2,2), mai=c(1.37,1.37,.7,.3))
#par(cex=1.25, mai=c(1.37,1.37,.7,.3) )
t <- 1464
proc_data$genes$gene_name[t]
x <- HTS_data$methyl_region[[t]][,1]
y <- HTS_data$methyl_region[[t]][,3]/HTS_data$methyl_region[[t]][,2]
xl <- seq(-1, 1, length = 250)
plot(x, y, col = "red3", pch = 21, ylim = c(0,1), xlim = c(-1,1), lwd = 2,
     xlab = NA, ylab = NA, cex.axis = 1.1,  xaxt = "n")
mtext(side = 1, "genomic region", line = 3, cex = 1.1)
mtext(side = 2, "methylation level", line = 3, cex = 1.1)
axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
title(main="Observations", line = 1, cex.main=1.3)
# lines(x = xl, y = eval_probit_function(out_mean$basis, xl, out_mean$W_opt[t,]), 
#       col = 'coral', lwd = 3, lty = 2)
# lines(x = xl, y = eval_probit_function(out_prof$basis, xl, out_prof$W_opt[t,]), 
#       col = 'red2', lwd = 3)

# dev.print(pdf, "../figures/plekhh3-2.pdf", width = 10.5, height = 7)
# dev.off()

cols <- c("firebrick", "cornflowerblue", "coral",  "darkolivegreen4", "#E69F00")
x.coord = seq(-1, 1, length = 250)

## RBF Centres
plot(NULL, xlim = c(-1,1), ylim = c(0,1), type = "n", xaxt = "n", xlab = NA, ylab = NA, cex.axis = 1.1)
mtext(side = 1, "genomic region", line = 3, cex = 1.1)
mtext(side = 2, "methylation level", line = 3, cex = 1.1)
axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
title(main="RBF centres", line = 1, cex.main=1.3)
points(basis_prof$mus, rep(0, basis_prof$M), col=cols[1:basis_prof$M], pch=19)
for (i in 1:basis_prof$M) {
  points(x.coord, pnorm(apply(X = as.matrix(x.coord), 
                        MARGIN = 1,
                        FUN = rbf_basis,
                        mus = basis_prof$mus[i],
                        gamma = basis_prof$gamma)), 
         type="l", lwd=2, col=cols[i])
}

## Fit RBFs using BPR
plot(NULL, xlim = c(-1,1), ylim = c(0,1), type = "n", xaxt = "n", xlab = NA, ylab = NA, cex.axis = 1.1)
mtext(side = 1, "genomic region", line = 3, cex = 1.1)
mtext(side = 2, "methylation level", line = 3, cex = 1.1)
axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
title(main="Learn coefficients using BPR model", line = 1, cex.main=1.3)
points(basis_prof$mus, rep(0, basis_prof$M), col=cols[1:basis_prof$M], pch=19)
for (i in 1:basis_prof$M) {
  evals <- pnorm(out_prof$W_opt[t,1] * apply(X = as.matrix(x.coord), 
                                             MARGIN = 1,
                                             FUN = rbf_basis,
                                             mus = basis_prof$mus[i],
                                             gamma = basis_prof$gamma) + 
                   out_prof$W_opt[t,i+1] * apply(X = as.matrix(x.coord), 
                                                 MARGIN = 1,
                                                 FUN = rbf_basis,
                                                 mus = basis_prof$mus[i],
                                                 gamma = basis_prof$gamma))
  points(x.coord, evals, type="l", lwd=2, col=cols[i])
}


## Linear combination
plot(x, y, col = "red3", pch = 21, ylim = c(0,1), xlim = c(-1,1), lwd = 2,
     xlab = NA, ylab = NA, cex.axis = 1.1,  xaxt = "n")
mtext(side = 1, "genomic region", line = 3, cex = 1.1)
mtext(side = 2, "methylation level", line = 3, cex = 1.1)
axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
title(main="Linear combination of RBFs", line = 1, cex.main=1.3)
lines(x = xl, y = eval_probit_function(out_prof$basis, xl, out_prof$W_opt[t,]), 
      col = 'red2', lwd = 3)

# #
# # DOCK5 methylation profile
# #
# which(proc_data$genes$gene_name == "DOCK5")
# xl <- seq(-1, 1, 0.001)
# 
# par(cex=1.25, mai=c(1.37,1.37,.7,.3) )
# t <- 4847 
# 
# t <- 1234
# x <- HTS_data$methyl_region[[t]][,1]
# y <- HTS_data$methyl_region[[t]][,3]/HTS_data$methyl_region[[t]][,2]
# xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
# plot(x, y, col = "red3", pch = 21, ylim = c(0,1), xlim = c(-1,1), lwd = 2,
#      xlab = NA, ylab = NA, cex.axis = 1.3)
# mtext(side = 1, "genomic region x", line = 3, cex = 1.8)
# mtext(side = 2, "methylation level", line = 3, cex = 1.8)
# title(main="Gene DOCK5", line = 1, cex.main=1.6)
# lines(x = xl, y = eval_probit_function(out_mean$basis, xl, out_mean$W_opt[t,]), 
#       col = 'coral', lwd = 3, lty = 2)
# lines(x = xl, y = eval_probit_function(out_prof$basis, xl, out_prof$W_opt[t,]), 
#       col = 'red2', lwd = 3)
