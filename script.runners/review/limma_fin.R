library(limma)
library(ggplot2)
library(RCurl)

getURL('https://raw.githubusercontent.com/pampernickel/subtyping/master/script.runners/routineFuncs.R', ssl.verifypeer = F) -> script
eval(parse(text=script))

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
sapply(cores.Iqbal, function(x) 
  which(colnames(sub.list.batch1$exprs) %in% x)) -> sub.list.batch1$groups

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
sapply(1:4, function(x) 
  which(sub.list$groups %in% x)) -> groups
setLimmaRun(sub.list$exprs.mat, groups) -> limma.res.t
save(limma.res.t, file="./r.data.files/results/gep/limma_tenomic.rda")

load("./r.data.files/external/sub.list.fin_genes.rda")
dropSamples(sub.list.batch1, which(sub.list.batch1$groups %in% 0)) -> sub.list.batch1
sapply(1:4, function(x) which(sub.list.batch1$groups %in% x)) -> g
setLimmaRun(sub.list.batch1$exprs, g) -> limma.res.i
save(limma.res.i, file="./r.data.files/results/gep/limma_iqbal.rda")

# ::: limma res comparisons
load("./r.data.files/results/gep/limma_tenomic.rda")
load("./r.data.files/results/gep/limma_iqbal.rda")
compareLimmaRes(limma.res.t, limma.res.i, "t") -> comp
plotComp(comp, "TENOMIC", "GSE58445")
compareLimmaRes(limma.res.t, limma.res.i, "logFC") -> comp
plotComp(comp, "TENOMIC", "GSE58445")

# ::: Check 1
# check group assignments of TENOMIC samples that were included in the Iqbal
# cohort; clear workspace of all variables as a precaution
rm(list = setdiff(ls(), lsf.str()))
load("./r.data.files/results/gep/limma_iqbal.rda")
ti <- read.csv("./annotations/tenomic_iqbal.csv")
gsub("TENOMIC", "CRE_TENOMIC_", ti$N..Biobase.TENOMIC..CRE_TENOMIC) -> ti
unlist(strsplit(ti, " "))[grep("TENOMIC", unlist(strsplit(ti, " ")))] -> tids

# ::: Check 2
# Iqbal samples in tenomic vs same background as before
rm(list = setdiff(ls(), lsf.str()))
load("./r.data.files/results/gep/limma_iqbal.rda")
ti <- read.csv("./annotations/tenomic_iqbal.csv")
gsub("TENOMIC", "CRE_TENOMIC_", ti$N..Biobase.TENOMIC..CRE_TENOMIC) -> ti
unlist(strsplit(ti, " "))[grep("TENOMIC", unlist(strsplit(ti, " ")))] -> tids

load("./r.data.files/second_proc/sub.list.fin_genes.rda")
dropSamples(sub.list, which(sub.list$groups %in% 0)) -> sub.list
sapply(unique(sub.list$groups), function(x) 
  intersect(which(sub.list$tenomic %in% tids), 
            which(sub.list$groups %in% x))) -> g
sapply(1:length(g), function(x) 
  setdiff(which(sub.list$groups %in% x), g[[x]])) -> e
setLimmaRun(sub.list$exprs.mat, g, e) -> limma.res.tiq
compareLimmaRes(limma.res.tiq, limma.res.i, "t") -> comp
plotComp(comp, "TENOMIC_IN_IQBAL", "GSE58445")

# ::: Check 3
# Quick check with possibility of Iqbal group label switch
# TBX21 and GATA3 marker expression level checks
load("./r.data.files/external/sub.list.fin_genes.rda")
dropSamples(sub.list.batch1, which(sub.list.batch1$groups %in% 0)) -> sub.list.batch1
sapply(unique(sub.list.batch1$groups), function(x) which(sub.list.batch1$groups %in% x)) -> g
plotGenes(sub.list.batch1$exprs, 
          g, c("GATA3", "LEF1", "CCR4", "TNFRSF8", "CXCL12"), 
          type=c("boxplot", "heatmap"))

load(file="./r.data.files/second_proc/sub.list.fin_genes.rda")
dropSamples(sub.list, which(sub.list$groups %in% 0)) -> sub.list
sapply(unique(sub.list$groups), function(x) 
  which(sub.list$groups %in% x)) -> g
colnames(sub.list$exprs.mat) <- sub.list$tenomic
plotGenes(sub.list$exprs.mat, 
          g, c("GATA3", "LEF1", "CCR4", "TNFRSF8", "CXCL12"), 
          type="boxplot")