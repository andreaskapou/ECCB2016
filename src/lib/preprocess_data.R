#' Function for preprocessing and filtering RNA-Seq data and keeping only the
#' corresponding promoter region for BS-Seq data.
#'
preprocess_data <-function(HTS_data,  is_fpkm = FALSE, max_outl = 600,
                           gene_expr_thresh = FALSE, gene_outl_thresh = FALSE,
                           gene_log2_transf = FALSE){

  # Get the FPKM values from RNA-Seq data, instead of default which is for
  # showing results in the genome browser.
  obs <- HTS_data$methyl_region
  if (is_fpkm){
    Y   <- as.numeric(HTS_data$rna_data$gene_fpkm)
    max_outl = 300
  }else{
    Y   <- HTS_data$rna_data$gene_expr
  }
  genes <- HTS_data$rna_data

  # Option to discard all unexpressed genes
  if (gene_expr_thresh){
    ind <- which(Y == 0)
    Y <- Y[-ind]
    obs <- obs[-ind]
    genes <- genes[-ind]
  }

  # Option to discard possible outliers / noisy data
  if (gene_outl_thresh){
    ind <- which(Y > max_outl)
    Y <- Y[-ind]
    obs <- obs[-ind]
    genes <- genes[-ind]
  }

  # Option to log-transform gene expression data
  if (gene_log2_transf){
    Y <- Y + 0.1
    Y <- log2(Y)
  }

  return(list(Y = Y, obs = obs, genes = genes))
}
