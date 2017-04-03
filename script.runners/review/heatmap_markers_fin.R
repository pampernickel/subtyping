library(gplots)

getURL('https://raw.githubusercontent.com/pampernickel/subtyping/master/script.runners/routineFuncs.R', ssl.verifypeer = F) -> script
eval(parse(text=script))

load("./r.data.files/second_proc/sub.list.fin_genes.rda") # limma_fin.r for processing
dropSamples(sub.list, which(sub.list$groups==0)) -> sub.list

# marker file
markers <- read.csv("./data/markers.csv")
markers[-which(duplicated(markers$marker)),] -> markers

# arrange df by group, and by markers
orderBy(sub.list, "groups") -> sub.list
sub.list$exprs.mat[which(rownames(sub.list$exprs.mat) %in% markers$marker),] -> df
binarizeMarkers(df, sub.list$groups, sub.list$tenomic, markers, 0.7)

load("./r.data.files/external/sub.list.fin_genes.rda") # limma_fin.r for processing
dropSamples(sub.list.batch1, which(sub.list.batch1$groups==0)) -> sub.list.batch1
orderBy(sub.list.batch1, "groups") -> sub.list.batch1
sub.list.batch1$exprs[which(rownames(sub.list.batch1$exprs) %in% markers$marker),] -> df
binarizeMarkers(df, sub.list.batch1$groups, colnames(sub.list.batch1$exprs), markers, 0.8)