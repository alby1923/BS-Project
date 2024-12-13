library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('ultimate_dataset.rds')

#DEPENDANCY GRAPH --------------
library(BDgraph)

target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
res_list <- list()

for (target in target_variables){
  load(paste0("Residuals_", target, ".Rdata"))
  res_list[[target]] <- res
}

res_data <- do.call(rbind, res_list)

res_data <- t(res_data)

bdgraph = bdgraph.mpl(res_data,iter = 20000)
summary(bdgraph)


