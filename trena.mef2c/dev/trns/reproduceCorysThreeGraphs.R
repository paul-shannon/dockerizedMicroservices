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

if(!exists("tbl.dhs"))
   load("../../trena/data/tbl.dhs.RData")

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
runTests <- function()
{
   test_makeGeneModel.footprints()
   test_makeGeneModel.encodeDHS()
   test_makeGeneModel.allDNA()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
makeGeneModel <- function(bindingSiteSource, tbl.roi, pfms=NA, pwmMinMatch=95, mtx, tfMappingSource,
                          targetGene, orderByColumn="pcaMax", solverInclusivenessCutoff=1.0)
{
   stopifnot(bindingSiteSource %in% c("footprints", "encodeDHS", "allDNA"))

   gr.regions <- GRanges(tbl.roi)

   if(bindingSiteSource == "footprints"){
      gr.fp <- GRanges(tbl.fp)
      tbl.overlaps <- as.data.frame(findOverlaps(gr.regions, gr.fp, type="any"))
      if(nrow(tbl.overlaps) == 0)
        return(list(trn=data.frame(), bindingSites=data.frame()))
      colnames(tbl.overlaps) <- c("roi", "fp")
      tbl.roiSites <- tbl.fp[tbl.overlaps$fp, ]
      if("geneSymbol" %in% colnames(tbl.roiSites)){
         column.to.remove <- grep("geneSymbol", colnames(tbl.roiSites))
         tbl.roiSites <- tbl.roiSites[, -column.to.remove]
         }
      } # footprints
   else if(bindingSiteSource == "encodeDHS"){
     gr.dhs <- GRanges(tbl.dhs)
     tbl.overlaps <- as.data.frame(findOverlaps(gr.regions, gr.dhs, type="any"))
     if(nrow(tbl.overlaps) == 0)
        return(list(trn=data.frame(), bindingSites=data.frame()))
     colnames(tbl.overlaps) <- c("roi", "dhs")
     tbl.roi.dhs <- tbl.dhs[tbl.overlaps$dhs, ]
     if(all(is.na(pfms)))
        pfms <- as.list(query(query(MotifDb("jaspar2018")), "hsapiens"))
      mm <- MotifMatcher(genomeName="hg38", pfms)
      tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.roi.dhs, pwmMatchMinimumAsPercentage=pwmMinMatch)
      if(nrow(tbl.motifs) == 0){
        printf("--- no match of pfms in supplied regions at %d%%", pwmMinMatch)
        return(list(trn=data.frame(), bindingSites=data.frame()))
        }
      tbl.roiSites <- tbl.motifs
      colnames(tbl.roiSites)[grep("motifStart", colnames(tbl.roiSites))] <- "start"
      colnames(tbl.roiSites)[grep("motifEnd", colnames(tbl.roiSites))] <- "end"
      tbl.roiSites$length <- 1 + tbl.roiSites$end - tbl.roiSites$start
      tbl.roiSites$distance.from.tss <- TSS - tbl.roiSites$start
      } # encodeDHS

   else if(bindingSiteSource == "allDNA"){
     if(all(is.na(pfms)))
        pfms <- as.list(query(query(MotifDb("jaspar2018")), "hsapiens"))
      mm <- MotifMatcher(genomeName="hg38", pfms)
      tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.roi, pwmMatchMinimumAsPercentage=pwmMinMatch)
      if(nrow(tbl.motifs) == 0){
        printf("--- no match of pfms in supplied regions at %d%%", pwmMinMatch)
        return(list(trn=data.frame(), bindingSites=data.frame()))
        }
      tbl.roiSites <- tbl.motifs
      colnames(tbl.roiSites)[grep("motifStart", colnames(tbl.roiSites))] <- "start"
      colnames(tbl.roiSites)[grep("motifEnd", colnames(tbl.roiSites))] <- "end"
      tbl.roiSites$length <- 1 + tbl.roiSites$end - tbl.roiSites$start
      tbl.roiSites$distance.from.tss <- TSS - tbl.roiSites$start
      } # allDNA


   tbl.roiSites$shortMotif <-  unlist(lapply(strsplit((tbl.roiSites$motifName), "-"),
                                             function(tokens) tokens[length(tokens)]))

   tbl.roiSites <- associateTranscriptionFactors(MotifDb, tbl.roiSites, source=tfMappingSource, expand.rows=TRUE)
   colnames(tbl.roiSites)[grep("geneSymbol", colnames(tbl.roiSites))] <- "tf"

   tfs <- unique(tbl.roiSites$tf)

   tfs.known <- intersect(tfs, rownames(mtx))   # 9 -> 7

   solverNames <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   solver <- EnsembleSolver(mtx, targetGene=targetGene,
                            candidateRegulators=tfs.known, solverNames, geneCutoff = solverInclusivenessCutoff)
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
test_makeGeneModel.footprints <- function()
{
   printf("--- test_makeGeneModel.footprints")

      # a small model
   tbl.roi <- data.frame(chrom="chr5", start=TSS, end=TSS+2000, stringsAsFactors=FALSE)

    #--------------------------------------------------------------------------------
    # use pre-calculated footprints, 2kb region, conservative mapping of motifs to tf
    # all tfs returned
    #--------------------------------------------------------------------------------

   x0 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi,
                       mtx=mtx.tcx,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="rfScore",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x0)), c("bindingSites", "trn"))
   tbl.trn.0 <- x0$trn
   tbl.bindingSites.0 <- x0$bindingSites
   tfs.0 <- tbl.trn.0$tf
   checkTrue(all(tfs.0 %in% tbl.bindingSites.0$tf))
   checkEquals(dim(tbl.trn.0), c(6, 10))
   checkEquals(dim(tbl.bindingSites.0), c(6, 14))
   checkEquals(sort(tfs.0), c("RARA::RXRA", "RUNX2", "SPI1", "TEAD1", "TEAD3", "ZNF740"))

    #--------------------------------------------------------------------------------
    # same as x0, but use TFclass motif-to-tf mapping
    #--------------------------------------------------------------------------------

     # same request, but using TFclass motif-to-tf mapping
   x1 <- makeGeneModel(bindingSiteSource="footprints",
                       tbl.roi,
                       mtx=mtx.tcx,
                       tfMappingSource="TFclass",
                       targetGene="MEF2C",
                       orderByColumn="rfScore",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x1)), c("bindingSites", "trn"))
   tbl.trn.1 <- x1$trn
   tbl.bindingSites.1 <- x1$bindingSites
   tfs.1 <- tbl.trn.1$tf
   checkTrue(length(setdiff(tfs.0, tfs.1)) > 0)   # TFClass finds more and different tfs than does MotifDb
   checkTrue(all(tfs.1 %in% tbl.bindingSites.1$tf))
   checkEquals(dim(tbl.trn.1), c(16, 10))
   checkEquals(dim(tbl.bindingSites.1), c(28, 14))
   checkTrue(all(tfs.1 %in% c("PRRX1", "ARX", "HOPX", "RXRA", "RARA", "TEAD1", "ZNF740", "RHOXF1", "VSX1", "OTX1",
                               "OTX2", "TEAD2", "TEAD3", "TEAD4", "UNCX", "SPI1")))

} # makeGeneModel.footprints
#------------------------------------------------------------------------------------------------------------------------
test_makeGeneModel.encodeDHS <- function()
{
   printf("--- test_makeGeneModel.encodeDHS")

      # a larger region producing a small model, since no dhs regions are near TSS
   tbl.roi <- data.frame(chrom="chr5", start=TSS, end=TSS+8000, stringsAsFactors=FALSE)

    #--------------------------------------------------------------------------------
    # now use encode DHS, jaspar2018 human, conservative tf mapping, all tfs
    #--------------------------------------------------------------------------------

   x2 <- makeGeneModel(bindingSiteSource="encodeDHS",
                       tbl.roi,
                       mtx=mtx.tcx,
                       pfms=as.list(query(query(MotifDb, "jaspar2018"), "hsapiens")),
                       pwmMinMatch=97,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="rfScore",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x2)), c("bindingSites", "trn"))
   tbl.trn.2 <- x2$trn
   tfs.2 <- tbl.trn.2$tf
   tbl.bindingSites.2 <- x2$bindingSites
   motifNames.2 <- tbl.bindingSites.2$shortMotif
   checkTrue(all(tfs.2 %in% tbl.bindingSites.2$tf))

   checkEquals(sort(unique(tbl.bindingSites.2$shortMotif)),
               c("MA0036.1", "MA0077.1", "MA0080.1", "MA0098.1", "MA0719.1", "MA1120.1", "MA1125.1"))

   checkEquals(dim(tbl.trn.2), c(7, 10))
   checkEquals(dim(tbl.bindingSites.2), c(16, 17))
   checkTrue(all(tfs.2 %in% c("RHOXF1", "SOX13", "ZNF384", "SOX9", "GATA2", "SPI1", "ETS1")))

    #--------------------------------------------------------------------------------
    # as x2, but with motif match at 95
    #--------------------------------------------------------------------------------

   x3 <- makeGeneModel(bindingSiteSource="encodeDHS",
                       tbl.roi,
                       mtx=mtx.tcx,
                       pfms=as.list(query(query(MotifDb, "jaspar2018"), "hsapiens")),
                       pwmMinMatch=95,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="rfScore",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x3)), c("bindingSites", "trn"))
   tbl.trn.3 <- x3$trn
   tbl.bindingSites.3 <- x3$bindingSites
   tfs.3 <- sort(unique(tbl.trn.3$tf))
   checkTrue(all(tfs.2 %in% tfs.3))
   checkTrue(length(tfs.3) > length(tfs.2))
   checkEquals(tfs.3, c("ETS1", "GATA2", "MEIS1", "RBPJ", "RHOXF1", "SOX13", "SOX15", "SOX9", "SPI1", "ZNF384"))

   motifNames.3 <- sort(unique(tbl.bindingSites.3$shortMotif))
   checkTrue(all(motifNames.2 %in% motifNames.3))
   checkEquals(motifNames.3, c("MA0036.1", "MA0077.1", "MA0080.1", "MA0098.1", "MA0498.2", "MA0719.1", "MA1116.1",
                               "MA1120.1", "MA1125.1", "MA1152.1"))

   checkTrue(all(tbl.bindingSites.3$motifRelativeScore >= 0.95))

   checkEquals(dim(tbl.trn.3), c(10, 10))
   checkEquals(dim(tbl.bindingSites.3), c(22, 17))
   checkEquals(sort(tbl.trn.3$tf), sort(unique(tbl.bindingSites.3$tf)))

} # test_makeGeneModel.encodeDHS
#------------------------------------------------------------------------------------------------------------------------
test_makeGeneModel.allDNA <- function()
{
   printf("--- test_makeGeneModel.allDNA")

   tbl.roi <- data.frame(chrom="chr5", start=TSS, end=TSS+2000, stringsAsFactors=FALSE)

    #--------------------------------------------------------------------------------
    # now use allDNA, jaspar2018 human, conservative tf mapping, all tfs
    #--------------------------------------------------------------------------------

   x5 <- makeGeneModel(bindingSiteSource="allDNA",
                       tbl.roi,
                       mtx=mtx.tcx,
                       pfms=as.list(query(query(MotifDb, "jaspar2018"), "hsapiens")),
                       pwmMinMatch=100,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="pcaMax",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x5)), c("bindingSites", "trn"))
   tbl.trn.5 <- x5$trn
   tfs.5 <- tbl.trn.5$tf
   tbl.bindingSites.5 <- x5$bindingSites
   motifNames.5 <- tbl.bindingSites.5$shortMotif
   checkTrue(all(tfs.5 %in% tbl.bindingSites.5$tf))

   checkEquals(sort(unique(tbl.bindingSites.5$shortMotif)),
               c("MA0056.1", "MA0080.1", "MA0098.1", "MA0124.1", "MA0157.1"))

   checkEquals(dim(tbl.trn.5), c(5, 10))
   checkEquals(dim(tbl.bindingSites.5), c(7, 17))
   checkTrue(all(tfs.5 %in% c("FOXO3", "ETS1", "NKX3-1", "SPI1", "MZF1")))

    #--------------------------------------------------------------------------------
    # as x5, but with motif match at 97
    #--------------------------------------------------------------------------------

   x6 <- makeGeneModel(bindingSiteSource="allDNA",
                       tbl.roi,
                       mtx=mtx.tcx,
                       pfms=as.list(query(query(MotifDb, "jaspar2018"), "hsapiens")),
                       pwmMinMatch=97,
                       tfMappingSource="MotifDb",
                       targetGene="MEF2C",
                       orderByColumn="rfScore",
                       solverInclusivenessCutoff=1.0)

   checkEquals(sort(names(x6)), c("bindingSites", "trn"))
   tbl.trn.6 <- x6$trn
   tbl.bindingSites.6 <- x6$bindingSites
   tfs.6 <- sort(unique(tbl.trn.6$tf))
   checkTrue(all(tfs.5 %in% tfs.6))
   checkEquals(tfs.6, c("ETS1", "FOXC1", "FOXG1", "FOXO3", "FOXO4", "FOXO6", "GATA2", "MZF1",
                        "NFATC2", "NFATC3", "NFIC", "NFIX", "NKX3-1", "PRRX1", "SP1", "SPI1",
                        "TEAD3", "UNCX", "ZEB1", "ZNF354C"))


   motifNames.6 <- sort(unique(tbl.bindingSites.6$shortMotif))
   checkTrue(all(motifNames.5 %in% motifNames.6))
   checkEquals(motifNames.6, c("MA0032.1", "MA0036.1", "MA0036.3", "MA0056.1", "MA0079.2", "MA0080.1",
                               "MA0098.1", "MA0103.2", "MA0103.3", "MA0124.1", "MA0130.1", "MA0152.1",
                               "MA0157.1", "MA0157.2", "MA0161.1", "MA0613.1", "MA0625.1", "MA0671.1",
                               "MA0716.1", "MA0721.1", "MA0808.1", "MA0848.1", "MA0849.1"))

   checkTrue(all(tbl.bindingSites.6$motifRelativeScore >= 0.97))

   checkEquals(dim(tbl.trn.6), c(20, 10))
   checkEquals(dim(tbl.bindingSites.6), c(44, 17))
   checkEquals(sort(tbl.trn.6$tf), sort(unique(tbl.bindingSites.6$tf)))

} # test_makeGeneModel.allDNA
#------------------------------------------------------------------------------------------------------------------------
