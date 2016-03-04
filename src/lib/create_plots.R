#' Plot scatterplot with correlation between predicted and measured gene
#' expression levels for test data.
#'
plot_scatt_regr_test <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.35, 8.02)
    xlim=c(-3.65, 6.72)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }

  # Compute correlation
  r <- round(stats::cor(output$test_pred, output$test$y), 2)
  # Compute RMSE
  rmse   <- round(output$test_errors$rmse, 2)

  plot(output$test_pred,
       output$test$y,
       ylab=NA,
       xlab=NA,
       ylim=ylim,
       xlim=xlim,
       cex.axis=1,
       col="#0000ff42",
       pch=16,
       cex=1.3)
  mtext(side = 1, "predicted expression (log2)", line = 2.2, cex = 1.1)
  mtext(side = 2, "measured expression (log2)", line = 2.2, cex = 1.1)
  title(main=main_lab, line = 1.1, cex.main=1.45)
  # Perform simple linear regression using lm
  my_lm = lm(output$test$y ~ output$test_pred)
  abline(stats::coef(my_lm)[1],
         stats::coef(my_lm)[2],
         col="red", lty=4, lwd=3)
  # Create legend inside the figure
  legend(-4.4, 8.8, #"topleft",
         legend = paste0("r = ", r, "\n", "RMSE = ", rmse),
         bty = "n", cex = 1.2)
}


#' ggplot scatterplot with correlation between predicted and measured gene
#' expression levels for test data.
#' 
ggplot_scatt_regr_test <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.35, 8.02)
    xlim=c(-3.85, 7.02)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }
  
  # Compute correlation
  r <- round(stats::cor(output$test_pred, output$test$y), 2)
  # Compute RMSE
  rmse   <- round(output$test_errors$rmse, 2)
  # Perform simple linear regression using lm
  my_lm = lm(output$test$y ~ output$test_pred)
  
  out_plot <- data.frame(pred = output$test_pred, 
                         meas = output$test$y)
  
  
  gg <- ggplot(out_plot, aes(x = pred, y = meas)) + 
    geom_point(pch = 16, col = "#0000ff56", cex = 3) + 
    theme_bw() +
    theme(axis.title.x = element_text(color="black", size=19),
          axis.title.y = element_text(color="black", size=19),
          plot.title = element_text(face="bold", color = "black", size=23),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.5)) + 
    labs(x = "predicted expression (log2)", y = "measured expression (log2)", 
         title=main_lab) + 
    scale_x_continuous(breaks = seq(-2, 8, by = 2), limits=xlim) + 
    scale_y_continuous(breaks = seq(-2, 8, by = 2), limits=ylim) + 
    geom_abline(intercept = stats::coef(my_lm)[1],
                slope = stats::coef(my_lm)[2],
                col="red", lty = 4, lwd = 1) +
    geom_text(data = data.frame(), 
              aes(-3.85, 7.85, label = paste0("r = ", r, "\n", "RMSE = ", rmse)),
              size = 6, colour = "black", 
              hjust = 0)
  
  return(gg)
}


#' Plot scatterplot with correlation between predicted and measured gene
#' expression levels for train data.
#'
plot_scatt_regr_train <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.2, 8.02)
    xlim=c(-3.45, 6.72)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }

  # Compute correlation
  r <- round(stats::cor(output$train_pred, output$train$y), 3)
  # Compute RMSE
  rmse   <- round(output$train_errors$rmse, 3)

  plot(output$train_pred,
       output$train$y,
       ylab="measured expression (log2)",
       xlab="predicted expression (log2)",
       main=main_lab,
       ylim=ylim,
       xlim=xlim,
       cex.main=1.6,
       cex.lab=1.4,
       cex.axis=1.2,
       col="#0000ff35",
       pch=16,
       cex=1.3)

  # Plot the best linear line
  my_lm = lm(output$train$y ~ output$train_pred)
  abline(stats::coef(my_lm)[1],
         stats::coef(my_lm)[2],
         col="red", lty=4, lwd=3)

  # Create legend inside the figure
  legend("topleft",
         legend = paste0("r = ", r, "\n", "RMSE = ", rmse),
         bty = "n", cex = 1.1)
}


#' Plot scatterplot with correlation across cell lines
#'
plot_scatt_across_cell_line <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.35, 8.02)
    xlim=c(-3.85, 6.62)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }
  
  # Compute correlation
  r <- round(stats::cor(output$test_pred, output$test$y), 2)
  # Compute RMSE
  rmse   <- round(output$test_errors$rmse, 2)
  
  plot(output$test_pred,
       output$test$y,
       ylab=NA,
       xlab=NA,
       ylim=ylim,
       xlim=xlim,
       cex.axis=1,
       col="#0000ff42",
       pch=16,
       cex=1.1)
  mtext(side = 1, "predicted expression (log2)", line = 2.2, cex = 1.1)
  mtext(side = 2, "measured expression (log2)", line = 2.2, cex = 1.1)
  title(main=main_lab, line = 1.1, cex.main=1.45, font = 2)
  # Perform simple linear regression using lm
  my_lm = lm(output$test$y ~ output$test_pred)
  abline(stats::coef(my_lm)[1],
         stats::coef(my_lm)[2],
         col="red", lty=4, lwd=3)
  # Create legend inside the figure
  legend(4.7, -1.5, #"topleft",
         legend = paste0("r = ", r, "\n", "RMSE = ", rmse),
         bty = "n", cex = 1.2)
}



#' Plot scatterplot with correlation across cell lines
#'
ggplot_scatt_across_cell_line <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.35, 8.02)
    xlim=c(-3.85, 6.62)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }
  
  # Compute correlation
  r <- round(stats::cor(output$test_pred, output$test$y), 2)
  # Compute RMSE
  rmse   <- round(output$test_errors$rmse, 2)
  # Perform simple linear regression using lm
  my_lm = lm(output$test$y ~ output$test_pred)
  
  out_plot <- data.frame(pred = output$test_pred, 
                         meas = output$test$y)
  
  
  gg <- ggplot(out_plot, aes(x = pred, y = meas)) + 
    geom_point(pch = 16, col="#0000ff56", cex = 1.6) + 
    theme_bw() +
    theme(axis.title.x = element_text(color="black", size=13.5),
          axis.title.y = element_text(color="black", size=13.5),
          plot.title = element_text(face="bold", color = "black", size=16.5),
          axis.text = element_text(size = 12),
          panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.5)) + 
    labs(x = "predicted expression (log2)", y = "measured expression (log2)", 
         title=main_lab) + 
    scale_x_continuous(breaks = seq(-4, 6, by = 2), limits=xlim) + 
    scale_y_continuous(breaks = seq(-4, 8, by = 2), limits=ylim) + 
    geom_abline(intercept = stats::coef(my_lm)[1],
                slope = stats::coef(my_lm)[2],
                col="red", lty = 4, lwd = 1) +
    geom_text(data = data.frame(), 
              aes(4.82, -2.35, label = paste0("r = ", r, "\n", "RMSE = ", rmse)),
              size = 3.7, colour = "black", 
              hjust = 0)
  
  return(gg)
}


plot_confusion_corr_matrix <- function(out_GM, out_GM_K562, out_GM_H1, 
                                       out_K562, out_K562_GM, out_K562_H1,
                                       out_H1, out_H1_GM, out_H1_K562, 
                                       title_lab = "Correlation"){
  
  # ----------------------------------------
  # Compute correlations from GM
  # ----------------------------------------
  cor_gm <- round(cor(out_GM$test_pred, out_GM$test$y), 3)
  cor_gm_k562 <- round(cor(out_GM_K562$test_pred, out_GM_K562$test$y), 3)
  cor_gm_h1 <- round(cor(out_GM_H1$test_pred, out_GM_H1$test$y), 3)
  
  # ----------------------------------------
  # Compute correlations from K562
  # ----------------------------------------
  cor_k562 <- round(cor(out_K562$test_pred, out_K562$test$y), 3)
  cor_k562_gm <- round(cor(out_K562_GM$test_pred, out_K562_GM$test$y), 3)
  cor_k562_h1 <- round(cor(out_K562_H1$test_pred, out_K562_H1$test$y), 3)
  
  # ----------------------------------------
  # Compute correlations from H1-hESC
  # ----------------------------------------
  cor_h1 <- round(cor(out_H1$test_pred, out_H1$test$y), 3)
  cor_h1_gm <- round(cor(out_H1_GM$test_pred, out_H1_GM$test$y), 3)
  cor_h1_k562 <- round(cor(out_H1_K562$test_pred, out_H1_K562$test$y), 3)

  # ---------------------------------------
  # Confusion correlation matrix 
  # ---------------------------------------
  df <- data.frame(C1 = c("K562","K562","K562","GM12878","GM12878",
                          "GM12878","H1-hESC","H1-hESC","H1-hESC"),
                   C2 = c("K562","GM12878","H1-hESC","GM12878","K562",
                          "H1-hESC","H1-hESC","K562","GM12878"), 
                   r = c(cor_k562, cor_k562_gm, cor_k562_h1, 
                         cor_gm, cor_gm_k562, cor_gm_h1, 
                         cor_h1, cor_h1_k562, cor_h1_gm),
                   stringsAsFactors = TRUE)
  
  confusion_corr_mat <- ggplot() +
    geom_tile(aes(x = C2, y = C1, fill = r), 
              data = df, color = "black", size = 0.1) +
    theme(axis.text = element_text(size = 12), 
          plot.title = element_text(face="plain", size = 16.5),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 7.5)) +
    geom_text(aes(x = C2, y = C1, label = sprintf("%.2f", r)),
              data = df, size = 4.5, colour = "black") +
    scale_fill_gradient2(low = "#006400", mid = "#f2f6c3", 
                         high = "#cd0000", limits = c(0, 1)) +
    labs(x = "",y = "") +
    ggtitle(title_lab)
  
  return(confusion_corr_mat)
}


plot_cluster_prof <- function(model, basis, add_clust = FALSE, main_lab = "Clustered methylation profiles"){
  xs <- seq(-1,1,len=2000) # create some values
  plot(x=xs, y=eval_probit_function(basis, xs, model$w[,1]), xlim=c(-1,1), ylim=c(0,1),
       type="l", col="salmon3", lwd=3, xlab="promoter region", ylab="methylation level",
       main=main_lab)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,2]),
        col="black", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,3]),
        col="blue", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,4]),
        col="red3", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,5]),
        col="darkgreen", lwd=3)

  if (add_clust){
    lines(x=xs, y=eval_probit_function(basis, xs, model$w[,6]),
          col="cyan", lwd=3)
    lines(x=xs, y=eval_probit_function(basis, xs, model$w[,7]),
          col="darkgoldenrod1", lwd=3)
  }
}


plot_cluster_box <- function(gene_expr, add_clust = FALSE, main_lab = "Gene expression levels"){
  col <- c("salmon3", "black", "blue", "red3", "darkgreen")
  if (add_clust){
    col <- c(col, "cyan", "darkgoldenrod1")
  }
  boxplot(gene_expr, col=col, notch=T, xlab="Cluster K",
          ylab="expression level", main=main_lab)
}