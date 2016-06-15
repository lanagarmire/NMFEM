
# generate data -----------------------------------------------------------

library(org.Mm.eg.db)
library(devtools)

set.seed(62742)

sub1 <- replicate(5, (seq(6)-60)*200 + rpois(6,12000))/100
sub2 <- replicate(5, (seq(6,1)-60)*200 + rpois(6,12000))/100
bg   <- replicate(5, rpois(6,12000))/1000 - 5
fpkms <- t(exp(cbind(sub1, bg, sub2))-1)
fpkms[fpkms < 0] <- 0

load('data/mppi.rda')

rownames(fpkms) <- sample(as.matrix(mppi), nrow(fpkms))
colnames(fpkms) <- paste0('sample_', seq(ncol(fpkms)))

fpkms <- fpkms[,sample(ncol(fpkms))]

toy <- fpkms

use_data(toy, overwrite=T)

# basic analyze and visualization of the dataset --------------------------

library(NMF)
library(ggplot2)
library(reshape2)
library(dplyr)

rwl <- log(fpkms+1)

ren <- nmf(rwl, rank=2, nrun=30, method="brunet", seed=12345, .options="p30v3")
consensusmap(ren)

basis(ren)
coef(ren)

ggdat <- melt(rwl) %>% tbl_df
ggplot() +
  geom_raster(aes(x=Var2, y=Var1, fill=value), ggdat)

