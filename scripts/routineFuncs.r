`%ni%` = Negate(`%in%`)

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
