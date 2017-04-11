`%ni%` = Negate(`%in%`)

dropSamples <- function(list, ind){
  # list is of the form
  # found in 20130902.full.list.genes.rda
  if (length(ind)==0){
    stop("0 length index. No changes made to your list.", call.=FALSE)
  }
  # --- for each list element, remove the indicated index
  for (i in 1:length(list)){
    if (is.matrix(list[[i]])){
      # --- drop columns
      list[[i]][,-ind] -> list[[i]]
    } else {
      list[[i]][-ind] -> list[[i]]
    }
  }
  
  return(list)
}

orderBy <- function(list, by=c("exprs.mat", "tenomic", "affy.no", "sorted", "hybrid.date", "tet2", "idh", "dnmt", "tfh.stat", "perc.tumor", "cell.type", "pathos", "centre","subgroup","class","hybrid.date.corrected")){
  vect <- list[[which(names(list) %in% by)]]
  #print(vect)
  for (i in 1:length(list)){
    if (is.matrix(list[[i]])){
      list[[i]][,order(vect)] -> list[[i]]
    } else {
      list[[i]][order(vect)] -> list[[i]]
    }
  }
  return(list)
}

runLimma <- function(df, lab){
  ff <- rep(0,ncol(df))
  ff[which(lab %in% 1)] <- 1
  design <- model.matrix(~ -1+factor(ff))
  colnames(design) <-c("GRP1", "GRP2") 
  fit <- lmFit(df, design) 
  contrast.matrix <- makeContrasts( 
    GRP.1.2 = GRP2 - GRP1,
    levels=design
  )
  
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  
  toptab.GRP = topTable(fit2, coef=c(1), number = dim(df)[1], adjust="BH")
  return(toptab.GRP)
}

setLimmaRun <- function(mat, g, e=""){
  # mat: exprs.mat
  # g: list containing columns of indices of mat; the length of g corresponds
  # with the number of subgroups; in the case of length(g) == 2, this results
  # in a pairwise comparison; in the case of length(g) > 2, the default mode
  # would be of the form x vs. rest; alternative mode to follow;
  # e: exclusion list; should be of the same length as g
  if (!is.loaded("limma")) require(limma)
  lapply(1:length(g), function(x){
    lab <- rep(0, ncol(mat))
    lab[g[[x]]] <- 1
    if (e %ni% "" && length(g) == length(e)){
      lab[e[[x]]] <- NA
    } else {
      stop("Exclusion list should have the same length as the group argument.")
    }
    mat[,which(lab %in% c(0,1))] -> mat
    lab[which(lab %in% c(0,1))] -> lab
    return(runLimma(mat, lab))
  }) -> limma.res
  return(limma.res)
}

compareLimmaRes <- function(l1, l2, nn){
  if (length(l1) != length(l2)) stop("Lists do not have equal lengths.")
  lapply(1:length(l1), function(x){
    if (nrow(l1[[x]]) != nrow(l2[[x]])) stop("Data frames do not have an equal number of rows.")
    l1[[x]][order(rownames(l1[[x]])),which(colnames(l1[[x]]) %in% nn)] -> l1x
    l2[[x]][order(rownames(l2[[x]])),which(colnames(l2[[x]]) %in% nn)] -> l2x
    cbind.data.frame(l1x,l2x,rep(paste("Subgroup" , x, "vs. REST", sep=" "), length(l1x))) -> t
    return(t)
  }) -> res
  do.call("rbind", res) -> res
  colnames(res) <- c("V1", "V2", "comparison")
  return(res)
}

plotComp <- function(res, xlab, ylab){
  # plot function for compareLimmaRes() results
  if (!is.loaded("ggplot2")) require(ggplot2)
  ggplot(res, aes(x=V1, y=V2))+
    geom_point(size=0.5)+facet_wrap(~comparison)+
    theme_bw()+xlab(xlab)+ylab(ylab)+
    theme(axis.title = element_text(size=18),
          axis.text=element_text(size=16),
          strip.text=element_text(size=17))
}

# ::: other useful functions for checks
plotGenes <- function(mat, g, gene, type){
  # plots gene expression across groups
  # mat: expression matrix
  # g: groups
  # gene: list of genes of interest
  if (!is.loaded("ggplot2")) require(ggplot2)
  if (!is.loaded("reshape2")) require(reshape2)
  if (!is.loaded("gplots")) require(gplots)
  if (!is.loaded("RColorBrewer")) require(RColorBrewer)
  
  if (type == "boxplot"){
    lapply(1:length(g), function(x){
      melt(mat[which(rownames(mat) %in% gene),g[[x]]]) -> t
      cbind.data.frame(t, rep(paste("Group",x,sep="_"), nrow(t))) -> df
      colnames(df) <- c(colnames(df)[1:(ncol(df)-1)], "group")
      return(df)
    }) -> res
    do.call("rbind", res) -> res
    ggplot(res, aes(x=group, y=value))+
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15, height=0.01, alpha=0.5)+
      facet_wrap(~Var1)+
      theme_bw()
  } else if (type == heatmap){
    mat[which(rownames(mat) %in% gene),] -> t
    my.colors <- colorRampPalette(colorRampPalette(brewer.pal(11,"RdBu")[-c(4,5,7,8)])(50))
    
    # create RowSide colors based on g
    heatmap.2(t(t), scale="col", col=my.colors,
              key=TRUE, keysize=0.5, trace="none", 
              sepwidth=0.01, #colsep=1:ncol(assembled), sepcolor='white',
              margins = c(10,10),
              cexRow=1,cexCol=1)
  }
}


# ::: create binarized heatmap
binarizeMarkers <- function(df, groups, sampleNames, markers, sdFactor){
  if (!is.loaded("gplots")) require(gplots)
  sapply(markers$marker, function(x) 
    which(rownames(df) %in% x)) -> rord
  df[rord,] -> df
  
  # currently hard-coded
  color.hash <- cbind(c(1:4), c("red", "cyan", "magenta", "blue"))
  colnames(color.hash) <- c("group", "color")
  col.colors <- rep(NA, ncol(df))
  for (i in 1:4){
    col.colors[which(groups %in% i)] <- color.hash[which(color.hash[,1] %in% i),2]
  }
  
  row.colors <- rep(NA, nrow(df))
  for (i in 1:4){
    row.colors[which(markers$group %in% i)] <- color.hash[which(color.hash[,1] %in% i),2]
  }
  
  # binarization; currently uses the median as the baseline limit
  # for binarization; the restriction can be adjusted with respect to the
  # sd in a vector
  factor <- sdFactor
  apply(df, 1, function(x)
    median(x)+factor*sd(x)) -> medx
  
  df.bin <- matrix(0, nrow=nrow(df), ncol=ncol(df))
  for (i in 1:nrow(df)){
    df.bin[i,which(df[i,] > medx[i])] <- 1
  }
  colnames(df.bin) <- sampleNames
  rownames(df.bin) <- rownames(df)
  
  # then within each group, order acc. to highest expression of markers
  fin.ord <- c()
  for (i in 1:4){
    df.bin[which(rownames(df.bin) %in% 
                   markers$marker[which(markers$group %in% i)]),
           which(groups %in% i)] -> x
    colSums(x) -> s
    colnames(x)[order(s)] -> n.ord
    c(fin.ord, n.ord) -> fin.ord
  }
  
  as.numeric(sapply(fin.ord, function(x) 
    which(colnames(df.bin) %in% x))) -> fin.ord
  df.bin[,fin.ord] -> df.bin
  
  heatmap.2(df.bin, dendrogram="none", Rowv=F, Colv=F,
            trace="none", col=c("grey95", "grey50"),
            ColSideColors=col.colors,
            RowSideColors=row.colors,
            labCol="", cexRow=1.5, margin=c(7,7),
            colsep=getChangePoint(as.numeric(as.factor(col.colors))), 
            rowsep=getChangePoint(as.numeric(as.factor(row.colors))), sepcolor = "grey12")
}

getChangePoint <- function(d) {
  p <- cumsum(rle(d)$lengths)
  p[-length(p)]
}

# ::: functions for enrichment analyses
runRomer <- function(df, labs, gene.sets){
  # --- params: MsigDB signature set, iset, sub.list, vector for generating design matrix
  ff <- rep(0,ncol(df))
  ff[which(labs %in% 1)] <- 1
  design <- cbind(Intercept=1,Group=labs)
  
  # index: list corresponding to the rows of the gene sets
  lapply(gene.sets, function(x) 
    ids2indices(x, rownames(df))) -> indices
  
  # select signatures of interest: 
  rr <- list()
  for (i in 1:length(gene.sets)){
    print(i)
    romer(y=df, index=indices[[i]], design=design, set.statistic = "mean", 
          nrot=50000) -> rr[[i]]
  }
  
  lapply(rr, function(x) 
    x[setdiff(c(which(x[,2] <= 0.05),
                which(x[,3] <= 0.05)), which(x[,4] <= 0.05)),]) -> rr.sub
  list(rr, rr.sub) -> res
  names(res) <- c("all", "sub")
  return(res)
}

createGeneSets <- function(dir){
  # specify directory from which to create gene sets; this function version has been
  # designed to handle MSigDB signatures
  list.files(dir, full.names = T)[grep(".txt", list.files(dir, full.names = T))] -> f
  lapply(f, function(x){
    readLines(x) -> ff
    lapply(ff, function(x) unlist(strsplit(x, "\t"), recursive=F)) -> ff
    lapply(ff, function(x) x[1]) -> names(ff)
    lapply(ff, function(x) x[-c(1:2)]) -> names(ff)
  }) -> gene.sets
  sapply(strsplit(f, "\\/"), function(x) x[length(x)]) -> names(gene.sets)
  return(gene.sets)
}