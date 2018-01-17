library(trena)
library(trenaViz)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("gene.models")){
   f <- "../../trena/data/mef2c.tf.5kb.RData"
   stopifnot(file.exists(f))
   mef2c.model.names <- load(f)    # mef2c.tcx, mef2c.cer, mef2c.ros
   printf("loading cory's 3 wg models: %s", paste(mef2c.model.names, collapse=", "))
   gene.models <- list()
   for(name in mef2c.model.names){  # 3 data.frames, 20 x 12
     gene.models[[name]] <- eval(parse(text=name))
     }
  } # read in gene.models
#------------------------------------------------------------------------------------------------------------------------
# can I really recreate cory's old models?  for instance, how comparable is the rosmap data?

#------------------------------------------------------------------------------------------------------------------------
add.dimer.expression <- function(mtx, dimers)
{
   for(dimer in dimers){
      genes <- strsplit(dimer, "::")[[1]]
      if(all(genes %in% rownames(mtx))){
        printf("adding dimer expression for %s", dimer)
        dimer.vec <- rep(1, ncol(mtx))
        for(gene in genes){
           mtx.sub <- mtx[gene,, drop=FALSE]
           #browser()
           dimer.vec <- dimer.vec * mtx.sub[1,]
           } # for gene
        mtx <- rbind(mtx, dimer.vec)
        rownames(mtx)[nrow(mtx)] <- dimer
        } # if all genes found
      } # for dimer

   invisible(mtx)

} # add.dimer.expression
#------------------------------------------------------------------------------------------------------------------------
test_add.dimer.expression <- function()
{
  dimer <- "PPARG::RXRA"
  mtx.dimers <- add.dimer.expression(mtx.cer, dimer)
  vec.1 <- mtx.cer["PPARG",]
  vec.2 <- mtx.cer["RXRA",]
  vec.1.2 <- vec.1 * vec.2
  checkEqualsNumeric(vec.1.2, mtx.dimers[dimer,])

  dimer <- "SMAD2::SMAD3::SMAD4"
  mtx.dimers <- add.dimer.expression(mtx.cer, dimer)
  vec.1 <- mtx.cer["SMAD2",]
  vec.2 <- mtx.cer["SMAD3",]
  vec.3 <- mtx.cer["SMAD4",]
  vec.1.2.3 <- vec.1 * vec.2 * vec.3
  checkEqualsNumeric(vec.1.2.3, mtx.dimers[dimer,])

  dimer <- "PPARG::bogus"
  mtx.dimers <- add.dimer.expression(mtx.cer, dimer)
  checkEquals(dim(mtx.dimers), dim(mtx.cer))

} # test_add.dimer.expression
#------------------------------------------------------------------------------------------------------------------------
# following Cory's suggestion, the expression matrices used here are from whovian, obtained (10 jan 2018):
#    ssh pshannon@whovian ls -l /local/Cory/Alzheimers/synapse.windsorized
#       -rw-r--r--. 1 cfunk pricelab       174 Nov 28 17:33 amp_ad_targets
#       drwxr-sr-x. 2 cfunk pricelab        10 Nov 28 14:31 gene_lists
#       -rw-r--r--. 1 cfunk pricelab  30642960 Jul 20 14:51 mayo.cer.RData
#       -rw-r--r--. 1 cfunk pricelab  30729909 Jul 20 14:49 mayo.tcx.RData
#       -rw-r--r--. 1 cfunk pricelab  68951488 Jul 27 13:31 rosmap.winsorized.RData
#       -rw-r--r--. 1 cfunk pricelab  81348050 Jul 20 14:40 Scaled_Winsorized_MayoRNAseq_CER.csv
#       -rw-r--r--. 1 cfunk pricelab  81531402 Jul 20 14:39 Scaled_Winsorized_MayoRNAseq_TCX.csv
#       -rw-r--r--. 1 cfunk pricelab 178570307 Jul 27 13:25 Scaled_Winsorized_ROSMAP_Expression.csv
# the three RData files are now in ../../trena/data.  note inconsistency in names.
if(!exists("tbl.fp")){
   load("../../trena/data/tbl.fp.chr5.88615025-89052115.4sources.noDups.RData")
   dimers <- unique(grep("::", tbl.fp$geneSymbol, value=TRUE))
   }


if(!exists("mtx.cer")){
   load("../../trena/data/mayo.cer.RData")
   mtx.cer <- as.matrix(tbl)               # 15160   263
   fivenum.string <- paste(round(fivenum(mtx.cer), 3), collapse=" ")
   printf("mtx.cer: %d x %d, %s", nrow(mtx.cer), ncol(mtx.cer), fivenum.string)
   mtx.cer <- add.dimer.expression(mtx.cer, dimers)
   }

if(!exists("mtx.tcx")){
   load("../../trena/data/mayo.tcx.RData")
   mtx.tcx <- as.matrix(tbl)                # 15160   264
   fivenum.string <- paste(round(fivenum(mtx.tcx), 3), collapse=" ")
   printf("mtx.tcx: %d x %d, %s", nrow(mtx.tcx), ncol(mtx.tcx), fivenum.string)
   mtx.tcx <- add.dimer.expression(mtx.tcx, dimers)
   }

if(!exists("mtx.ros")){
   load("../../trena/data/rosmap.winsorized.RData")
   mtx.ros <- as.matrix(tbl)                # 14228   632
   fivenum.string <- paste(round(fivenum(mtx.ros), 3), collapse=" ")
   printf("mtx.ros: %d x %d, %s", nrow(mtx.ros), ncol(mtx.ros), fivenum.string)
   mtx.ros <- add.dimer.expression(mtx.ros, dimers)
   }

if(!exists("tbl.snp"))
    load("../../trena/data/tbl.snp.hg38.score-ref-alt.RData")

#------------------------------------------------------------------------------------------------------------------------
if(!exists("trena"))
   trena <- Trena("hg38")

if(!exists("tv")){
   tv <- trenaViz(8000:8100)
   setGenome(tv, "hg38")
   }

#------------------------------------------------------------------------------------------------------------------------
targetGene <- "MEF2C"
TSS <- 88904257
fp.roi <- list(chrom="chr5", start=88615025, end=89052115)  # "chr5:88,615,025-89,052,115"
#roi.
roi.5kb  <- sprintf("chr5:%d-%d", TSS-5000, TSS+5000)
foi.eqtl <- sprintf("chr5:%d-%d", min(tbl.snp$pos)-50000, max(tbl.snp$pos)+50000)
#------------------------------------------------------------------------------------------------------------------------
makeGeneModel <- function(bindingSiteSource, tbl.roi, motifs.provided=FALSE, mtx, tfMappingSource,
                          targetGene, orderByColumn="pcaMax")
{
   stopifnot(bindingSiteSource %in% c("footprints")) #, "encodeDHS", "allDNA"))

   gr.regions <- GRanges(tbl.roi)

   if(bindingSiteSource == "footprints"){
      gr.fp <- GRanges(tbl.fp)
      tbl.overlaps <- as.data.frame(findOverlaps(gr.regions, gr.fp))
      colnames(tbl.overlaps) <- c("roi", "fp")
      tbl.roiSites <- tbl.fp[tbl.overlaps$fp, ]
      if("geneSymbol" %in% colnames(tbl.roiSites)){
         column.to.remove <- grep("geneSymbol", colnames(tbl.roiSites))
         tbl.roiSites <- tbl.roiSites[, -column.to.remove]
         }
      } # footprints

      # add in TFClass tfs
   #colnames(tbl.roiSites)[grep("geneSymbol", colnames(tbl.roiSites))] <- "tf"
   tbl.roiSites$shortMotif <-  unlist(lapply(strsplit((tbl.roiSites$motifName), "-"), "[", 4))

   tbl.roiSites <- associateTranscriptionFactors(MotifDb, tbl.roiSites, source=tfMappingSource, expand.rows=TRUE)
   colnames(tbl.roiSites)[grep("geneSymbol", colnames(tbl.roiSites))] <- "tf"

   tfs <- unique(tbl.roiSites$tf)

   tfs.known <- intersect(tfs, rownames(mtx))   # 9 -> 7

   solverNames <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   solver <- EnsembleSolver(mtx, targetGene=targetGene,
                            candidateRegulators=tfs.known, solverNames, geneCutoff = 0.5)
   tbl.model <- run(solver)
   colnames(tbl.model)[grep("gene", colnames(tbl.model))] <- "tf"
   tbl.bindingSitesInModel <- subset(tbl.roiSites, tf %in% tbl.model$tf)

   coi <- c("chrom", "start", "end", "motifName", "strand", "length", "distance.from.tss", "tf")
   dups <- which(duplicated(tbl.bindingSitesInModel[, coi]))
   if(length(dups) > 0)
      tbl.bindingSitesInModel <- tbl.bindingSitesInModel[-dups,]

   tbl.bindingSiteCounts <- as.data.frame(table(tbl.bindingSitesInModel$tf))
   colnames(tbl.bindingSiteCounts) <- c("tf", "bindingSites")

   tbl.model <- merge(tbl.model, tbl.bindingSiteCounts, by="tf")
   tbl.out <- tbl.model[order(tbl.model[, orderByColumn], decreasing = TRUE),]
   rownames(tbl.out) <- NULL
   list(trn=tbl.out, bindingSites=tbl.bindingSitesInModel)

} # makeGeneModel
#------------------------------------------------------------------------------------------------------------------------
test_makeGeneModel <- function()
{
   printf("--- test_makeGeneModel")

      # a small model
   tbl.roi <- data.frame(chrom="chr5", start=TSS, end=TSS+2000, stringsAsFactors=FALSE)

   x0 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi,
                       motifs.provided=TRUE,
                       mtx=mtx.tcx,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="rfScore")

   checkEquals(sort(names(x0)), c("bindingSites", "trn"))
   tbl.trn <- x0$trn
   tbl.bindingSites <- x0$bindingSites
   checkEquals(dim(tbl.trn), c(6, 10))
   checkEquals(dim(tbl.bindingSites), c(6, 14))
   checkEquals(sort(tbl.trn$tf), sort(unique(tbl.bindingSites$tf)))
   checkTrue(all(tbl.trn$tf %in%  c("TEAD1", "ZNF740", "TEAD3", "RARA::RXRA", "RUNX2", "SPI1")))

     # same request, but using TFclass motif-to-tf mapping
   x1 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi,
                       motifs.provided=TRUE,
                       mtx=mtx.tcx,
                       tfMappingSource="TFclass",
                       targetGene="MEF2C",
                       orderByColumn="rfScore")

   checkEquals(sort(names(x1)), c("bindingSites", "trn"))
   tbl.trn <- x1$trn
   tbl.bindingSites <- x1$bindingSites
   checkEquals(dim(tbl.trn), c(8, 10))
   checkEquals(dim(tbl.bindingSites), c(13, 14))
   checkEquals(sort(tbl.trn$tf), sort(unique(tbl.bindingSites$tf)))
   checkTrue(all(tbl.trn$tf %in%  c("ARX", "PRRX1", "HOPX", "RARA", "RXRA", "ZNF740", "TEAD1", "TEAD2")))


} # test_makeGeneModel
#------------------------------------------------------------------------------------------------------------------------
