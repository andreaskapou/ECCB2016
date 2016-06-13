if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(ggplot2)
library(cowplot)
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)

res <- "../files/cor_rbf_regr_models_MonMar071856.RData"
load(res)
print(model_name)
round(gm_out_prof[[1]]$test_errors$pcc, 2)
round(gm_out_prof[[2]]$test_errors$pcc, 2)
round(gm_out_prof[[3]]$test_errors$pcc, 2)
round(gm_out_prof[[4]]$test_errors$pcc, 2)

print(model_name)
round(k562_out_prof[[1]]$test_errors$pcc, 2)
round(k562_out_prof[[2]]$test_errors$pcc, 2)
round(k562_out_prof[[3]]$test_errors$pcc, 2)
round(k562_out_prof[[4]]$test_errors$pcc, 2)

print(model_name)
round(h1_out_prof[[1]]$test_errors$pcc, 2)
round(h1_out_prof[[2]]$test_errors$pcc, 2)
round(h1_out_prof[[3]]$test_errors$pcc, 2)
round(h1_out_prof[[4]]$test_errors$pcc, 2)
