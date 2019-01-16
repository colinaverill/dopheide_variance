#Calculate intra-core variance!!!
rm(list=ls())
source('paths.r')
library(DirichletReg)

#load data.----
d <- readRDS(variance_analysis_product.path)
map <- d$map
tax <- d$tax.reads

#subset to 1.5g direct extraction.
map <- map[map$DNA_extraction_method == '1.5 g direct',]
#Choose one sample.
map <- map[map$Sample == 'P_Psoil',]

#Work up phylum data to relative abundance.----
#1. subset to phyla in at least 50% of samples.
#2. put everything else in other, convert to relative abundances.
#3. Do re-sampling to generate variance curve.
test <- tax$phylum
check <- list()
for(i in 1:ncol(test)){
  check[[i]] <- sum(test[,i] > 0) / nrow(test)
}
check <- unlist(check)
names(check) <- colnames(test)
test <- test[, colnames(test) %in% names(check[check > 0.5])]

#deal with zeros, convert to relative abundance, subset to samples in map.----
test <- data.frame(test)
test$other <- d$seq_depth - rowSums(test)
test <- test+1
test <- test/rowSums(test)
test <- test[rownames(test) %in% map$run,]

#CALCULATE A VARIANCE!
test$Y <- DR_data(test[,1:ncol(test)])
mod <- DirichReg(Y ~ 1, data = test)
mult <- sum(mod$fitted.values$alpha[1,])
