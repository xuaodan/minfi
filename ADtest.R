
library(FlowSorted.Blood.450k)
#referenceset
mset <- preprocessQuantile(FlowSorted.Blood.450k)
#get user samples
#AD=read.table('C:\\Users\\XIAO\\Desktop\\minfi\\data\\AD.beta')
#AD=read.table('C:\\Users\\XIAO\\Desktop\\minfi\\data\\AD0.beta')
#AD=read.table('C:\\Users\\XIAO\\Desktop\\minfi\\data\\normal.beta')
AD=read.table('C:\\Users\\XIAO\\Desktop\\minfi\\data\\benign.beta')
#AD=read.table('C:\\Users\\XIAO\\Desktop\\minfi\\data\\malignant.beta')
ADset=as.matrix(AD)
ADset=makeGenomicRatioSetFromMatrix(ADset)

#cpg match
referenceSet=mset[rownames(ADset)]
# defination of split
splitit <- function(x) {
  split(seq(along=x), x)
}

p <- getBeta(referenceSet)
pd <- as.data.frame(colData(referenceSet))
#extract celltypes we need
keep <- which(pd$CellType %in% c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
pd <- pd[keep,]
p <- p[,keep]
#make cell type a factor
pd$CellType <- factor(pd$CellType, levels = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
#f-test between different cell types
ffComp <- rowFtests(p, pd$CellType)
#get rowmeans value of different cell types
prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
r <- matrixStats::rowRanges(p)
#generate all evaluation into one table/ f-test,rowmeans,max,min
compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
names(compTable)[1] <- "Fstat"
names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
# t-test for cell types
tIndexes <- splitit(pd$CellType)
tstatList <- lapply(tIndexes, function(i) {
  x <- rep(0,ncol(p))
  x[i] <- 1
  return(rowttests(p, factor(x)))
})
#get probeList
numProbes=50
probeList <- lapply(tstatList, function(x) {
  y <- x[x[,"p.value"] < 1e-8,]
  yUp <- y[order(y[,"dm"], decreasing=TRUE),]
  yDown <- y[order(y[,"dm"], decreasing=FALSE),]
  c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
})

trainingProbes <- unique(unlist(probeList))
trainingProbes=trainingProbes[!is.na(trainingProbes)]
p <- p[trainingProbes,]

#
pMeans <- colMeans(p)
names(pMeans) <- pd$CellType
form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
#call functions
tmp <- get("validationCellType", env=environment(fun=estimateCellCounts))(Y = p, pheno = phenoDF, modelFix = form)
coefEsts <- tmp$coefEsts
counts2 <- get("projectCellType", env=environment(fun=estimateCellCounts))(getBeta(ADset)[rownames(coefEsts),], coefEsts) 

