library(trena)
library(trenaViz)
print(load("data/mtx.protectedAndExposedLowVarianceGenesEliminated.RData"))
trena <- Trena("hg38")
source.1 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_16"
source.2 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_20"
source.3 <- "postgres://bddsrds.globusgenomics.org/skin_hint_16"
source.4 <- "postgres://bddsrds.globusgenomics.org/skin_hint_20"
sources <- c(source.1, source.2, source.3, source.4)
names(sources) <- c("well_16", "well_20", "hint_16", "hint_20")
sources <- sources[4]
tss <- 50201632

#------------------------------------------------------------------------------------------------------------------------
buildModelForRegion <- function(roi)
{
   cls <- parseChromLocString(roi)   # a trena function
   x <- getRegulatoryChromosomalRegions(trena, cls$chrom, cls$start, cls$end, sources, targetGene, tss)
   names(x) <- names(sources)
   x2 <- lapply(names(x), function(name) {tbl <-x[[name]]; if(nrow(tbl) >0) tbl$db <- name; return(tbl)})
   tbl.reg <- do.call(rbind, x2)
   rownames(tbl.reg) <- NULL
   tbl.reg <- unique(tbl.reg[grep("Hsapiens-jaspar2016", tbl.reg$motifName, ignore.case=TRUE),])
   tbl.reg <- tbl.reg[order(tbl.reg$motifStart),]
   tbl.reg$fpName <- paste(tbl.reg$motifName, tbl.reg$db, sep="_")

   solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   #solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
   colnames(tbl.reg)[c(2,3)] <- c("start", "end")
   tbl.reg.tfs <- associateTranscriptionFactors(MotifDb, tbl.reg, source="MotifDb", expand.rows=FALSE)
   tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.reg.tfs, mtx)
   tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]
   tbl.reg.tfs <- subset(tbl.reg.tfs, geneSymbol %in% tbl.geneModel$gene)
   list(model=tbl.geneModel, regions=tbl.reg.tfs)

} # buildModelForRegion
#------------------------------------------------------------------------------------------------------------------------
# tss: 50201632
loc_1 <- 'chr17:50,201,610-50,201,700'
loc_2 <- 'chr17:50,201,720-50,202,700'
targetGene <- "COL1A1"

left <- buildModelForRegion(loc_1)
right<- buildModelForRegion(loc_2)

model <- list(left=left, right=right)
g <- buildMultiModelGraph(targetGene, model)

xCoordinate.span <- 1500
g.lo <- addGeneModelLayout(g, xPos.span=xCoordinate.span)
g.json <- trenaViz:::.graphToJSON(g.lo)

PORT.RANGE <- 8000:8020
if(!exists("tv")){
   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
   setGenome(tv, "hg38")
   }

setGraph(tv, g.lo, names(model))
setStyle(tv, "../hostDir/trenaVizStyle.js")
