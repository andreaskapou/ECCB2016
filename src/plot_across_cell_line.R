# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(cowplot)
library(ggplot2)
library(earth)
library(randomForest)
library(e1071)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)


# -----------------------------------------
# Initialize variables
# -----------------------------------------
filename <- "../files/across_cell_line_corr_ThuMar031511.RData"

# -----------------------------------------
# Load saved data
# -----------------------------------------
load(filename)


# -----------------------------------------
# Create input features for each cell line
# -----------------------------------------

# GM features
GM_prof_W <- data.frame(x = GM_out_prof$W_opt,
                        y = GM_proc_data$Y)
GM_mean_W <- data.frame(x = GM_out_mean$W_opt,
                        y = GM_proc_data$Y)
# K562 features
K562_prof_W <- data.frame(x = K562_out_prof$W_opt,
                          y = K562_proc_data$Y)
K562_mean_W <- data.frame(x = K562_out_mean$W_opt,
                          y = K562_proc_data$Y)
# H1 features
H1_prof_W <- data.frame(x = H1_out_prof$W_opt,
                        y = H1_proc_data$Y)
H1_mean_W <- data.frame(x = H1_out_mean$W_opt,
                        y = H1_proc_data$Y)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# ----------------------------------------
# Predict from GM to all other cell lines
# ----------------------------------------
message("Predicting from GM ...")
GM_prof_pred_K562 <- predict_model_gex(model      = GM_out_prof$gex_model,
                                       test       = K562_prof_W,
                                       is_summary = TRUE)
GM_mean_pred_K562 <- predict_model_gex(model      = GM_out_mean$gex_model,
                                       test       = K562_mean_W,
                                       is_summary = TRUE)

GM_prof_pred_H1 <- predict_model_gex(model      = GM_out_prof$gex_model,
                                     test       = H1_prof_W,
                                     is_summary = TRUE)
GM_mean_pred_H1 <- predict_model_gex(model      = GM_out_mean$gex_model,
                                     test       = H1_mean_W,
                                     is_summary = TRUE)


# ----------------------------------------
# Predict from K562 to all other cell lines
# ----------------------------------------
message("Predicting from K562 ...")
K562_prof_pred_GM <- predict_model_gex(model      = K562_out_prof$gex_model,
                                       test       = GM_prof_W,
                                       is_summary = TRUE)
K562_mean_pred_GM <- predict_model_gex(model      = K562_out_mean$gex_model,
                                       test       = GM_mean_W,
                                       is_summary = TRUE)

K562_prof_pred_H1 <- predict_model_gex(model      = K562_out_prof$gex_model,
                                       test       = H1_prof_W,
                                       is_summary = TRUE)
K562_mean_pred_H1 <- predict_model_gex(model      = K562_out_mean$gex_model,
                                       test       = H1_mean_W,
                                       is_summary = TRUE)


# ----------------------------------------
# Predict from H1 to all other cell lines
# ----------------------------------------
message("Predicting from H1-hESC ...")
H1_prof_pred_GM <- predict_model_gex(model      = H1_out_prof$gex_model,
                                     test       = GM_prof_W,
                                     is_summary = TRUE)
H1_mean_pred_GM <- predict_model_gex(model      = H1_out_mean$gex_model,
                                     test       = GM_mean_W,
                                     is_summary = TRUE)

H1_prof_pred_K562 <- predict_model_gex(model      = H1_out_prof$gex_model,
                                       test       = K562_prof_W,
                                       is_summary = TRUE)
H1_mean_pred_K562 <- predict_model_gex(model      = H1_out_mean$gex_model,
                                       test       = K562_mean_W,
                                       is_summary = TRUE)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# ----------------------------------------
# GM prediction output
# ----------------------------------------
out_GM_prof_K562 <- list(test_pred = GM_prof_pred_K562$test_pred,
                         test = list(y = K562_proc_data$Y))
out_GM_mean_K562 <- list(test_pred = GM_mean_pred_K562$test_pred,
                         test = list(y = K562_proc_data$Y))

out_GM_prof_H1 <- list(test_pred = GM_prof_pred_H1$test_pred,
                       test = list(y = H1_proc_data$Y))
out_GM_mean_H1 <- list(test_pred = GM_mean_pred_H1$test_pred,
                       test = list(y = H1_proc_data$Y))


# ----------------------------------------
# K562 prediction output
# ----------------------------------------
out_K562_prof_GM <- list(test_pred = K562_prof_pred_GM$test_pred,
                         test = list(y = GM_proc_data$Y))
out_K562_mean_GM <- list(test_pred = K562_mean_pred_GM$test_pred,
                         test = list(y = GM_proc_data$Y))

out_K562_prof_H1 <- list(test_pred = K562_prof_pred_H1$test_pred,
                         test = list(y = H1_proc_data$Y))
out_K562_mean_H1 <- list(test_pred = K562_mean_pred_H1$test_pred,
                         test = list(y = H1_proc_data$Y))


# ----------------------------------------
# H1-hESC prediction output
# ----------------------------------------
out_H1_prof_GM <- list(test_pred = H1_prof_pred_GM$test_pred,
                       test = list(y = GM_proc_data$Y))
out_H1_mean_GM <- list(test_pred = H1_mean_pred_GM$test_pred,
                       test = list(y = GM_proc_data$Y))

out_H1_prof_K562 <- list(test_pred = H1_prof_pred_K562$test_pred,
                         test = list(y = K562_proc_data$Y))
out_H1_mean_K562 <- list(test_pred = H1_mean_pred_K562$test_pred,
                         test = list(y = K562_proc_data$Y))


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# ---------------------------------------
# Make plots
# ---------------------------------------
out_GM_prof_K562$test_errors <- GM_prof_pred_K562$test_errors
out_GM_mean_K562$test_errors <- GM_mean_pred_K562$test_errors
corr_prof <- ggplot_scatt_across_cell_line(output = out_GM_prof_K562, 
                                           main_lab = expression(Profile~GM12878 %->% K562), 
                                           is_margins = TRUE)

corr_mean <- ggplot_scatt_across_cell_line(output = out_GM_mean_K562, 
                                           main_lab = expression(Mean~GM12878 %->% K562), 
                                           is_margins = TRUE)

k562_file <- "../files/corr_K562_WedMar021931.RData"
gm_file   <- "../files/corr_GM_WedMar021943.RData"
h1_file   <- "../files/corr_H1_WedMar021830.RData"

load(k562_file)
K562_profile <- out_prof[[1]]
K562_mean <- out_mean[[1]]

load(gm_file)
GM_profile <- out_prof[[1]]
GM_mean <- out_mean[[1]]

load(h1_file)
H1_profile <- out_prof[[1]]
H1_mean <- out_mean[[1]]


confus_prof <- plot_confusion_corr_matrix(GM_out_prof, out_GM_prof_K562, out_GM_prof_H1, 
                                          K562_out_prof, out_K562_prof_GM, out_K562_prof_H1,
                                          H1_out_prof, out_H1_prof_GM, out_H1_prof_K562, 
                                          K562_profile, GM_profile, H1_profile,
                                          title_lab = "Methylation Profiles Correlation")


confus_mean <- plot_confusion_corr_matrix(GM_out_mean, out_GM_mean_K562, out_GM_mean_H1, 
                                          K562_out_mean, out_K562_mean_GM, out_K562_mean_H1,
                                          H1_out_mean, out_H1_mean_GM, out_H1_mean_K562, 
                                          K562_mean, GM_mean, H1_mean,
                                          title_lab = "Mean Methylation Correlation")


across_cell_plot <- plot_grid(confus_prof, corr_prof, confus_mean,  corr_mean, labels = c("A", "C", "B", ""), 
          label_size = 25, ncol = 2, nrow = 2)

# save_plot("../figures/corr-across-cell-lines.pdf", across_cell_plot, ncol = 2, nrow = 2,
#           base_aspect_ratio = 1.35)