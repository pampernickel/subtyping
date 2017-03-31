library(limma)
library(ggplot2)

source('./scripts/routineFuncs.r')

#::: TENOMIC feature selection (from version mapped to symbols)
load("./r.data.files/second_proc/sub.list.normalized.adjusted.combat.cov_QC.noSex.mapped.symbol.entrezID.w.groups.Rdata") #161 samples
cbind(sub.list$tenomic, sub.list$groups) -> group

# for each unique gene id, keep the one with the most variability based on sd, rudimentary check
sapply(unique(rownames(sub.list$exprs.mat)), function(x) 
  length(which(rownames(sub.list$exprs.mat) %in% x))) -> check

sapply(unique(rownames(sub.list$exprs.mat)), function(x){
  sub.list$exprs.mat[which(rownames(sub.list$exprs.mat) %in% x),] -> loc
  if (is.data.frame(loc) | is.matrix(loc)){
    loc[which(apply(loc, 1, function(x) sd(x)) %in% max(apply(loc, 1, function(x) sd(x))))[1],] -> loc
  }
  return(loc)
}) -> mat
t(mat) -> mat
mat -> sub.list$exprs.mat
# save(sub.list, file="./r.data.files/second_proc/sub.list.fin_genes.rda")

#::: Iqbal feature selection (from version mapped to symbols)
load("./r.data.files/external/sub.list.batch1.Rdata")
load("./r.data.files/external/cores.Iqbal.Rdata")
sub.list.batch1$groups <- rep(0, ncol(sub.list.batch1$exprs))
for (i in 1:length(cores.Iqbal)){
  sub.list.batch1$groups[which(colnames(sub.list.batch1$exprs) %in% cores.Iqbal[[i]])] <- i
}

# map to gene symbol and perform gene selection
load("./r.data.files/annots/ann_with_symbol.rda")
sapply(rownames(sub.list.batch1$exprs), function(x){
  return(ann_with_symbol$hgnc_symbol[which(ann_with_symbol$probeset_id %in% x)])
}) -> gn
rownames(sub.list.batch1$exprs) <- gn

sapply(unique(rownames(sub.list.batch1$exprs)), function(x){
  sub.list.batch1$exprs[which(rownames(sub.list.batch1$exprs) %in% x),] -> loc
  if (is.data.frame(loc) | is.matrix(loc)){
    loc[which(apply(loc, 1, function(x) sd(x)) %in% max(apply(loc, 1, function(x) sd(x))))[1],] -> loc
  }
  return(loc)
}) -> mat
t(mat) -> mat

# select same genes as in TENOMIC
mat[which(rownames(mat) %in% rownames(sub.list$exprs.mat)),] -> mat
mat -> sub.list.batch1$exprs
# save(sub.list.batch1, file="./r.data.files/external/sub.list.fin_genes.rda")

# ::: TENOMIC and Iqbal diff expr analyses
# stick to core samples, perform x vs REST comparison
load("./r.data.files/second_proc/sub.list.fin_genes.rda")
dropSamples(sub.list, which(sub.list$groups %in% c(0, NA))) -> sub.list
sapply(unique(sub.list$groups), function(x) 
  which(sub.list$groups %in% x)) -> groups
setLimmaRun(sub.list$exprs.mat, groups) -> limma.res.t
# save(limma.res.t, file="./r.data.files/results/gep/limma_tenomic.rda")

load("./r.data.files/external/sub.list.fin_genes.rda")
dropSamples(sub.list.batch1, which(sub.list.batch1$groups %in% 0)) -> sub.list.batch1
sapply(unique(sub.list.batch1$groups), function(x) which(sub.list.batch1$groups %in% x)) -> g
setLimmaRun(sub.list.batch1$exprs, g) -> limma.res.i
# save(limma.res.i, file="./r.data.files/results/gep/limma_iqbal.rda")

# ::: limma res comparisons
load("./r.data.files/results/gep/limma_tenomic.rda")
load("./r.data.files/results/gep/limma_iqbal.rda")
compareLimmaRes(limma.res.t, limma.res.i, "t") -> comp
plotComp(comp, "TENOMIC", "GSE58445")
compareLimmaRes(limma.res.t, limma.res.i, "logFC") -> comp
plotComp(comp, "TENOMIC", "GSE58445")

# ::: Check 1
# check group assignments of TENOMIC samples that were included in the Iqbal
# cohort
rm(list=ls())
source('./scripts/routineFuncs.r')
load("./r.data.files/results/gep/limma_iqbal.rda")
ti <- read.csv("./annotations/tenomic_iqbal.csv")
gsub("TENOMIC", "CRE_TENOMIC_", ti$N..Biobase.TENOMIC..CRE_TENOMIC) -> ti
unlist(strsplit(ti, " "))[grep("TENOMIC", unlist(strsplit(ti, " ")))] -> tids

# ::: Check 2
# Iqbal samples in tenomic vs same background as before
rm(list=ls())
source('./scripts/routineFuncs.r')
load("./r.data.files/results/gep/limma_iqbal.rda")
ti <- read.csv("./annotations/tenomic_iqbal.csv")
gsub("TENOMIC", "CRE_TENOMIC_", ti$N..Biobase.TENOMIC..CRE_TENOMIC) -> ti
unlist(strsplit(ti, " "))[grep("TENOMIC", unlist(strsplit(ti, " ")))] -> tids

load("./r.data.files/second_proc/sub.list.fin_genes.rda")
sub.list$groups[which(sub.list$tenomic %in% tids)] -> g
sub.list$tenomic[which(sub.list$tenomic %in% tids)] -> t

limma.res.tiq <- list()
for (i in 1:4){
  lab <- rep(0, length(sub.list$tenomic))
  lab[which(sub.list$groups %in% i)] <- NA
  lab[which(sub.list$tenomic %in% 
              t[which(g %in% i)])] <- 1
  sub.list$lab <- lab
  dropSamples(sub.list, which(is.na(sub.list$lab))) -> loc
  runLimma(loc$exprs.mat, loc$lab) -> limma.res.tiq[[i]]
}

df.sum <- matrix(0, nrow=0, ncol=3)
colnames(df.sum) <- c("GSE58445", "TENOMIC", "comparison")
for (i in 1:length(limma.res.i)){
  limma.res.i[[i]][order(rownames(limma.res.i[[i]])),] -> iq
  limma.res.tiq[[i]][order(rownames(limma.res.tiq[[i]])),] -> ten
  cbind(iq$t, ten$t, rep(paste("Subgroup" ,i, "vs. REST", sep=" "), 
                         nrow(ten))) -> t
  colnames(t) <- colnames(df.sum)
  rbind(df.sum, t) -> df.sum
}

as.data.frame(df.sum) -> df.sum
for (i in 1:2){
  as.numeric(as.character(df.sum[,i])) -> df.sum[,i]
}
ggplot(df.sum, aes(x=TENOMIC, y=GSE58445))+
  geom_point(size=0.5)+facet_wrap(~comparison)+
  theme_bw()

# ::: Check 3
# check -- by subsampling *within* each group -- which samples are
# responsible for the shape