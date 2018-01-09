library(rzmq)
library(jsonlite)
library(trena)
library(trenaViz) # used here only for buildMultiModelGraph and addGeneLayout
library(RUnit)
PORT <- 5548
targetGene <- "MEF2C"
targetGene.tss <- 88904257
#------------------------------------------------------------------------------------------------------------------------
# load expression matrices on startup
#
expression.matrix.files <- list(rosmapFrontalCortex="../data/rosmap.fcx.RData",
                                mayoTemporalCorex="../data/mayo.tcx.RData")
expression.matrices <- list()

# stash all computed objects here, to save return trips from jupyter and, especially, awkward JSON conversions of data.frams
cache <- new.env(parent=emptyenv())

for(matrix.name in names(expression.matrix.files)){
    file.name <- expression.matrix.files[[matrix.name]]
    printf("about to load expression matrix %s from %s", matrix.name, file.name)
    mtx.name <- load(file.name)
    eval(parse(text=sprintf("mtx <- %s", mtx.name)))
    expression.matrices[[matrix.name]] <- mtx
    }

# three regulatory models provided by cory using rosmap? and +/- 5kb promoter
mef2c.model.names <- load("../data/mef2c.tf.5kb.RData")
printf("loading cory's 3 wg models: %s", paste(mef2c.model.names, collapse=", "))
gene.models <- list()
for(name in mef2c.model.names){
   gene.models[[name]] <- eval(parse(text=name))
   }

# eqtl variants using hg38 coordinates
printf("loading variant data: %s", load("../data/tbl.snp.hg38.score-ref-alt.RData"))

# load up previously obtained footprints for a generously large region, larger than the span of
# the mef2c eqtl snps.  uses all 4 sources: hint|wellington, seed 20/16
printf("loading footprints: %s", load("../data/tbl.fp.chr5.88615025-89052115.4sources.noDups.RData"))
# we will check incoming fp(roi) requests against what we already have
fp.roi <- list(chrom="chr5", start=88615025, end=89052115)  # "chr5:88,615,025-89,052,115"

#------------------------------------------------------------------------------------------------------------------------
# precalculate two dhs/encode tables:
#
#  1) tbl.dhs:  all the open chromatin regions in fp.roi
#  2) tbl.dhsMotif:  all the jaspar2016/hsapiens motifs above threshold in those dhs regions
#
# 1) tbl.dhs <- getEncodeDHSRegions("hg38", "wgEncodeRegDnaseClustered", "chr5", fp.roi$start, fp.roi$end)
#    colnames(tbl.dhs) <- c("chrom", "start", "end", "count", "score")
#    save(tbl.dhs, file="../data/tbl.dhs.RData")
#
# 2)
#    sources <- c("encodeHumanDHS")
#    trena <- Trena("hg38")
#    tbl.dhsMotifs <- getRegulatoryChromosomalRegions(trena, fp.roi$chrom, fp.roi$start, fp.roi$end, sources,
#                                                     targetGene, targetGene.tss)
#    tbl.dhsMotifs <- tbl.dhsMotifs[[1]]
#    dups <- which(duplicated(tbl.dhsMotifs[, c(1,2,3,5)]))
#     if(length(dups) > 0)
#        tbl.dhsMotifs <- tbl.dhsMotifs[-dups,]
#     tbl.dhsMotifs <- associateTranscriptionFactors(MotifDb, tbl.dhsMotifs, source="MotifDb", expand.rows=TRUE)
#     dim(tbl.dhsMotifs) # 28872 9
#     save(tbl.dhsMotifs, file="../data/tbl.dhsMotifs.RData")
#------------------------------------------------------------------------------------------------------------------------
# load the open chromatin (dhs) table
printf("loading dhs regions: %s", load("../data/tbl.dhs.RData"))
printf("loading motifs in dhs regions: %s", load("../data/tbl.dhsMotifs.RData"))
#------------------------------------------------------------------------------------------------------------------------
# obtain all motifs, independent of chromatin state and footprints, across the entire extended region
#    pfms <- as.list(query(query(MotifDb, "jaspar2016"), "hsapiens"))
#    mm <- MotifMatcher(genomeName="hg38", pfms)
#    tbl.allDNAMotifs <- findMatchesByChromosomalRegion(mm, as.data.frame(fp.roi), pwmMatchMinimumAsPercentage=85)
#    tbl.allDNAMotifs <- associateTranscriptionFactors(MotifDb, tbl.allDNAMotifs, source="MotifDb", expand.rows=TRUE)
# pull out the standard columns
#    tbl.allDNAMotifs <- tbl.allDNAMotifs[, c("chrom", "motifStart", "motifEnd", "motifName", "strand", "motifRelativeScore", "geneSymbol")]
# rename start and end
#    colnames(tbl.allDNAMotifs)[2:3] <- c("start", "end")
#    save(tbl.allDNAMotifs, file="../data/tbl.allDNAMotifs.RData")   # 233600 x 14
printf("loading motifs found across all DNA: %s", load("../data/tbl.allDNAMotifs.RData"))
#------------------------------------------------------------------------------------------------------------------------
# load in the normalized tcx matrix whose samples (columns) correspond to pheno and geno information.
# prepared like this:
#     ~/s/work/priceLab/AD/expression.matrix.prep
#     prepareMatrices.R starts off like this:
#        load("ampADMayo.64253genes.278samples.RData")
#        stopifnot(dim(mtx) == c(64253, 278))
#        covariates.file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
#     and produces
#        print(load("~/s/work/priceLab/AD/expression.matrix.prep/prepped.tcx.matrices.RData"))
#           "mtx.tcx"     "mtx.tcx.ctl" "mtx.tcx.ad"
#	fivenum(mtx.tcx); fivenum(mtx.tcx.ctl); fivenum(mtx.tcx.ad)
#          [1] 0.000000e+00 2.269386e+00 1.384377e+01 4.185981e+01 6.446054e+05
#          [1] -9.965784  1.152255  3.758768  5.350882 19.298057
#          [1] -9.965784  1.233155  3.814262  5.387233 17.440484
#       dim(mtx.tcx); dim(mtx.tcx.ctl); dim(mtx.tcx.ad)
#          [1] 18281   278
#          [1] 18281    80
#          [1] 18281    84
#       mtx.tcx.normalized <- asinh(mtx.tcx)
#       fivenum(mtx.tcx.normalized)
#          [1]  0.000000  1.558003  3.322285  4.427616 14.069541
#       head(colnames(mtx.tcx.normalized))
#          [1] "S11344_TCX" "S11316_TCX" "S11431_TCX" "S11341_TCX" "S11289_TCX" "S11327_TCX"
#       samples.mtx <- sub("^S", "", colnames(mtx.tcx.normalized))
#       samples.mtx <- sub("_TCX", "", samples.mtx)
#       colnames(mtx.tcx.normalized) <- samples.mtx
#
#       length(samples.mtx)  # [1] 278
#       length(samples.geno) # [1] 349
#       length(intersect(samples.mtx, samples.geno)) # [1] 263
#
#    covariates.file <- "~/s/work/priceLab/AD/expression.matrix.prep/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
#    tbl.covariates <- read.table(covariates.file, sep=",", as.is=TRUE, header=TRUE)
#    dim(tbl.covariates) # [1] 278   8
#    as.data.frame(table(tbl.covariates$Diagnosis))
#                     Var1 Freq
#       1               AD   84
#       2          Control   80
#       3 Pathologic Aging   30
#       4              PSP   84
#    tbl.covariates$sample <- gsub("_TCX", "", tbl.covariates$ID)
#    length(intersect(samples.geno, tbl.covariates$sample))  # [1] 263
#    tbl.pheno <- subset(tbl.covariates, sample %in% samples.geno)  # 263 9
#    as.data.frame(table(tbl.pheno$Diagnosis))
#                  Var1 Freq
#    1               AD   79
#    2          Control   73
#    3 Pathologic Aging   30
#    4              PSP   81
#    samples.ad <- subset(tbl.pheno, Diagnosis=="AD")$sample
#    samples.ctl <- subset(tbl.pheno, Diagnosis=="Control")$sample
#    cor(mtx.tcx.normalized["MEF2C", samples.ctl], mtx.tcx.normalized["NFATC3", samples.ctl])  # [1] -0.8153529
#    cor(mtx.tcx.normalized["MEF2C", samples.ad], mtx.tcx.normalized["NFATC3", samples.ad])    # [1] -0.7929382
#
#    save(tbl.pheno, mtx.geno, samples.ad, samples.ctl, mtx.tcx.normalized, file="~/github/dockerizedMicroservices/trena.mef2c/trena/data/mtx.tcx.pheno.geno.RData")
#    load("~/github/dockerizedMicroservices/trena.mef2c/trena/data/mtx.tcx.pheno.geno.RData")
#
#
printf("loading expression, pheno and geno: %s", paste(load("../data/mtx.tcx.pheno.geno.RData"), collapse=", "))
#------------------------------------------------------------------------------------------------------------------------
# interpret a data.frame as a list of data (a bare data.frame), columnames, and rownames
# this makes for easy reconstruction into a pandas dataframe in python.
dataFrameToPandasFriendlyList <- function(tbl)
{
   rownames = rownames(tbl)
   colnames = colnames(tbl)
   rownames(tbl) <- NULL

   list(rownames=rownames, colnames=colnames, tbl=tbl)

} # dataFrameToPandasFriendlyList
#------------------------------------------------------------------------------------------------------------------------
processWellStructuredMessage <- function(msg)
{
   if(msg$cmd == "summarizeExpressionMatrices"){
      tbl.summary.as.list <- summarizeExpressionMatrices()
      response <- list(cmd=msg$callback, status="success", callback="", payload=tbl.summary.as.list)
      }
   else if(msg$cmd == "getSessionInfo"){
      info <- as.character(sessionInfo())
      response <- list(cmd=msg$callback, status="success", callback="", payload=info)
      }
   else if(msg$cmd == "getExpressionMatrixNames"){
      info <- sort(names(expression.matrix.files))
      response <- list(cmd=msg$callback, status="success", callback="", payload=info)
      }
   else if(msg$cmd == "getModelNames"){
      info <- sort(mef2c.model.names)
      response <- list(cmd=msg$callback, status="success", callback="", payload=info)
      }
   else if(msg$cmd == "getModel"){
      modelName <- msg$payload
      key <- as.character(as.numeric(Sys.time()) * 100000)
      tbl <- gene.models[[modelName]]
      cache[[key]] <- tbl
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "getVariants"){
      thresholdScore <- msg$payload$minScore
      tbl.snp$filteringScore <- -log10(tbl.snp$CER_P)
      tbl.var <- subset(tbl.snp, filteringScore >= thresholdScore)[, c("chrom", "pos", "pos", "rsid", "filteringScore")]
      colnames(tbl.var) <- c("chrom", "start", "end", "id", "score")
      rownames(tbl.var) <- NULL
      key <- as.character(as.numeric(Sys.time()) * 100000)
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "getFootprintsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roi <- msg$payload$roi
      tbl.reg <- getFootprints(roi)
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.reg
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "getDHSRegionsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roi = msg$payload$roi
      tbl.dhs.roi <- getDHSRegions(roi)
      if(nrow(tbl.dhs.roi) == 0){
         response <- list(cmd=msg$callback, status="failure", callback="", payload="no overlapping regions");
         }
      else{
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.dhs.roi
         tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.dhs.roi)
         payload <- list(tbl=tbl.fp.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         } # else: some overlap found
      } # getDHsRegionsInRegion
   else if(msg$cmd == "getDHSMotifsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roi <- msg$payload$roi
      tbl.reg <- getDHSMotifs(roi)
      if(nrow(tbl.reg) == 0){
         response <- list(cmd=msg$callback, status="failure", callback="", payload="no motifs in region");
         }
      else{
         key <- as.character(as.numeric(Sys.time()) * 100000)
         cache[[key]] <- tbl.reg
         tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.reg)
         payload <- list(tbl=tbl.fp.as.list, key=key)
         response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
         }
      } # getDHSMotifsInRegion
   else if(msg$cmd == "findVariantsInModel"){
      stopifnot("modelName" %in% names(msg$payload))
      stopifnot("shoulder" %in% names(msg$payload))
      modelName <- msg$payload$modelName
      shoulder <- msg$payload$shoulder
      tbl.var <- findVariantsInModel(modelName, shoulder, "footprints")
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.var
      tbl.fp.as.list <- dataFrameToPandasFriendlyList(tbl.var)
      payload <- list(tbl=tbl.fp.as.list, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "listSharedData"){
      filenames <- list.files("/home/trena/sharedData")
      response <- list(cmd=msg$callback, status="success", callback="", payload=filenames);
      }
   else if(msg$cmd == "createTaggedDataFrame"){
      rows <- msg$payload$rows
      cols <- msg$payload$cols
      tbl <- as.data.frame(matrix(runif(32), nrow=rows, ncol=cols, dimnames=list(LETTERS[1:rows],  letters[1:cols])))
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl
      tbl.for.pandas <- dataFrameToPandasFriendlyList(tbl)
      payload <- list(tbl=tbl.for.pandas, key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "identifyTaggedDataFrame"){
      key <- msg$payload
      found <- key %in% names(cache)
      tbl.retrieved <- cache[[key]]
      printf("retrieved tbl with dimensions %d, %d", nrow(tbl.retrieved), ncol(tbl.retrieved))
      print(tbl.retrieved)
      sum <- sum(as.matrix(tbl.retrieved))
      payload <- list(found=found, sum=sum)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload);
      }
   else if(msg$cmd == "createGeneModel"){
      trena <- Trena("hg38")   # probably should create a single global instead
      payload <- msg$payload
      solver.names <- payload$solverNames
      printf("---- createGeneModel about to extract tbl.motifs from cache")
      key <- payload$tblRegulatoryRegionsCacheKey
      printf("    key: %s", key)
      found <- key %in% names(cache)
      printf("    found? %s", found)
      tbl.motifs <- cache[[key]]
      printf("    tbl.motifs: %d, %d", nrow(tbl.motifs), ncol(tbl.motifs))
      tfMap <- payload$tfMap
      stopifnot(all(c("targetGene", "matrixName") %in% names(msg$payload)))
      targetGene <- toupper(msg$payload$targetGene)
      matrixName <- msg$payload$matrixName
      #tbl.motifs <- read.table("/home/trena/sharedData/tbl.bed", sep="\t", as.is=TRUE, stringsAsFactors=FALSE)
      tbl.motifs <- tbl.motifs[, 1:5]
      colnames(tbl.motifs) <- c("chrom", "start", "end", "motifName")
      print(head(tbl.motifs))
      tbl.motifs$motifName <- sub("_well_16", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_well_20", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_hint_16", "", tbl.motifs$motifName)
      tbl.motifs$motifName <- sub("_hint_20", "", tbl.motifs$motifName)
      print(head(tbl.motifs))

      tbl.motifs.tfs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source=tfMap, expand.rows=FALSE)
      #print(head(tbl.motifs.tfs))

      mtx <- expression.matrices[[matrixName]]
      stopifnot(targetGene %in% rownames(mtx))
      print(dim(mtx))
      print(mean(mtx[targetGene,]))
      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs.tfs, mtx)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]
      key <- as.character(as.numeric(Sys.time()) * 100000)
      cache[[key]] <- tbl.geneModel
      print(tbl.geneModel)
      payload=list(tbl=dataFrameToPandasFriendlyList(tbl.geneModel), key=key)
      response <- list(cmd=msg$callback, status="success", callback="", payload=payload)
      } # createGeneModel
   else if(msg$cmd == "buildMultiModelGraph"){
      printf("----- buildMultiModelGraph")
      payload = msg$payload
      print(payload)
      targetGene <- payload$targetGene
      models <- payload$models
      models.decached <- list()
      for(name in names(models)){
         model <- models[[name]]
         tbl.geneModel <-  cache[[model$model]]
         tbl.regions <- cache[[model$regions]]
         tbl.regions <- subset(tbl.regions, geneSymbol %in% tbl.geneModel$gene)
         models.decached[[name]]$model <- tbl.geneModel
         models.decached[[name]]$regions <- tbl.regions
         } # for model
      printf("---- after decaching")
      print(models.decached)
      g <- trenaViz::buildMultiModelGraph(targetGene, models.decached)
      xCoordinate.span <- 1500
      g.lo <- trenaViz::addGeneModelLayout(g, xPos.span=xCoordinate.span)
      printf("---- g.lo")
      print(g.lo)
      g.json <- trenaViz:::.graphToJSON(g.lo)
      response <- list(cmd=msg$callback, status="success", callback="", payload=g.json)
      } # buildMultiModelGraph
   else{
      response <- list(cmd=msg$callback, status="success", callback="", payload="well-structured (unparsed) message")
      }

   response

} # processWellStructuredMessage
#------------------------------------------------------------------------------------------------------------------------
summarizeExpressionMatrices <- function()
{
   matrix.count <- length(expression.matrices)
   empty.vec <- rep(0, matrix.count)
   tbl.mtx <- data.frame(name=names(expression.matrices),
                         nrow=empty.vec,
                         ncol=empty.vec,
                         min=empty.vec,
                         q1=empty.vec,
                         median=empty.vec,
                         q3=empty.vec,
                         max=empty.vec,
                         stringsAsFactors=FALSE);
   rownames(tbl.mtx) <- tbl.mtx$name
   tbl.mtx <- tbl.mtx[, -1]

   for(name in rownames(tbl.mtx)){
      mtx <- expression.matrices[[name]]
      summary.stats <- fivenum(mtx)
      dimensions <- dim(mtx)
      tbl.mtx[name, c("min", "q1", "median", "q3", "max")] <- summary.stats
      tbl.mtx[name, c("nrow", "ncol")] <- dimensions
      } # for mtx.name

     # this result is destined for JSON and a python pandas dataframe
     # structure the data.frame as a 3-part list for easy uptake with pandas

   dataFrameToPandasFriendlyList(tbl.mtx)

} # summarizeExpressionMatrices
#------------------------------------------------------------------------------------------------------------------------
getFootprints <- function(roiString)
{
   roi <- parseChromLocString(roiString)   # a trena function
   if(roi$chrom == fp.roi$chrom &
      roi$start >= fp.roi$start &
      roi$end   <= fp.roi$end) {
      tbl.roi <- subset(tbl.fp, chrom==roi$chrom & start >= roi$start & end <= roi$end)
      return(tbl.roi)
      }

   trena <- Trena("hg38")
   source.1 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_16"
   source.2 <- "postgres://bddsrds.globusgenomics.org/brain_wellington_20"
   source.3 <- "postgres://bddsrds.globusgenomics.org/brain_hint_16"
   source.4 <- "postgres://bddsrds.globusgenomics.org/brain_hint_20"
   sources <- c(source.1, source.2, source.3, source.4)
   names(sources) <- c("well_16", "well_20", "hint_16", "hint_20")

   #sources <- sources[4];


   x <- getRegulatoryChromosomalRegions(trena, roi$chrom, roi$start, roi$end, sources, targetGene, targetGene.tss)
   print(1)
   names(x) <- names(sources)
   print(2)

      # append a column to each non-empty table, giving it the source name
   x2 <- lapply(names(x), function(name) {tbl <-x[[name]]; if(nrow(tbl) >0) tbl$db <- name; return(tbl)})
   print(3)

   tbl.reg <- do.call(rbind, x2)
   print(4)
   rownames(tbl.reg) <- NULL
   print(5)
      # be strict for now: just the 2016 jaspar human motifs
   tbl.reg <- unique(tbl.reg[grep("Hsapiens-jaspar2016", tbl.reg$motifName, ignore.case=TRUE),])
   tbl.reg <- tbl.reg[order(tbl.reg$motifStart),]
   printf("---- getFootprints, before associateTranscriptionFactors")
   print(head(tbl.reg))
   tbl.reg <- associateTranscriptionFactors(MotifDb, tbl.reg, source="MotifDb", expand.rows=TRUE)

   #printf("---- getFootprints, after associateTranscriptionFactors")
   #print(head(tbl.reg))

      # nasty hack, working around igv.js, npm, ipywidgets, docker, ???: last line in tbl is dropped
   colnames(tbl.reg)[2:3] <- c("start", "end")
   #printf("--- writing tbl.bed")
   #write.table(tbl.bed, file="/home/trena/sharedData/tbl.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
   #printf("--- after write")
   print(6)
   tbl.reg

} # getFootprints
#------------------------------------------------------------------------------------------------------------------------
getDHSRegions <- function(roiString)
{
   roi <- parseChromLocString(roiString)   # a trena function
   tbl.roi <- subset(tbl.dhs, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   return(tbl.roi)

} # getDHSRegions
#------------------------------------------------------------------------------------------------------------------------
getDHSMotifs <- function(roiString)
{
   roi <- parseChromLocString(roiString)   # a trena function
   tbl.roi <- subset(tbl.dhsMotifs, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   return(tbl.roi)

} # getDHSMotifs
#------------------------------------------------------------------------------------------------------------------------
findVariantsInModel <- function(modelName, shoulder, motifSources)
{
      #--------------------------------------------------------------------------------
      # clean up the snps dataframe, putting it in good bed file style
      # we then look, below, for intersections of these locations with
      # 3 sources for motifs:  bdds footprints, in encode dhs clustered, in all DNA
      #--------------------------------------------------------------------------------

   stopifnot(all(motifSources %in% c("footprints", "dhs", "allDNA")))

   tbl.snp.clean <- tbl.snp[, c("chrom", "pos", "pos")]
   colnames(tbl.snp.clean) <- c("chrom", "start", "end")
   gr.snp <- GRanges(tbl.snp.clean)

      # fill into these data.frames below as appropriate: if implied by motifSources parameter,
      # if any hits are found

   tbl.fpHits <- data.frame()
   tbl.dhsMotifHits <- data.frame()
   tbl.allMotifHits <- data.frame()

      #--------------------------------------------------------------------------------
      # do the bdds footprints firest
      #--------------------------------------------------------------------------------

   if("footprints" %in% motifSources){
      gr.fp <- GRanges(tbl.fp)
      tbl.fp.ov <- as.data.frame(findOverlaps(gr.snp, gr.fp, maxgap=shoulder))
      if(nrow(tbl.fp.ov) > 0){
         colnames(tbl.fp.ov) <- c("snp", "fp")
         tbl.fpHits <- cbind(tbl.snp[tbl.fp.ov$snp,], tbl.fp[tbl.fp.ov$fp, "geneSymbol", drop=FALSE ])
         tbl.fpHits <- unique(subset(tbl.fpHits, geneSymbol %in% gene.models[[modelName]]$gene))
         if(nrow(tbl.fpHits) > 0){
            tbl.fpHits$modelName <- modelName
            tbl.fpHits$shoulder <- shoulder
            tbl.fpHits$source <- "hint+wellington, 16+20 footprints"
            } # if motif's tf is in the gene model
         } # if motifs overlapped with footprints
      } # if footprints requested

      #--------------------------------------------------------------------------------
      # now the motifs found in encode open chromatin
      #--------------------------------------------------------------------------------

   if("dhs" %in% motifSources){
      gr.dhsMotifs <- GRanges(tbl.dhsMotifs)
      tbl.dhs.ov <- as.data.frame(findOverlaps(gr.snp, gr.dhsMotifs, maxgap=shoulder))
      if(nrow(tbl.dhs.ov) > 0){
         colnames(tbl.dhs.ov) <- c("snp", "dhsMotif")
         tbl.dhsMotifHits <- cbind(tbl.snp[tbl.dhs.ov$snp,],
                                   tbl.dhsMotifs[tbl.dhs.ov$dhsMotif, c("geneSymbol"), drop=FALSE])
         tbl.dhsMotifHits <- unique(subset(tbl.dhsMotifHits, geneSymbol %in% gene.models[[modelName]]$gene))
         if(nrow(tbl.dhsMotifHits) > 0){
            tbl.dhsMotifHits$modelName <- modelName
            tbl.dhsMotifHits$shoulder <- shoulder
            tbl.dhsMotifHits$source <- "encode DHS motifs"
            } # if hits in gene model
          } # if overlap
       } # if dhs requested

      #--------------------------------------------------------------------------------
      # now the motifs found in *all* dna
      #--------------------------------------------------------------------------------
   if("allDNA" %in% motifSources){
      gr.allDNAMotifs <- GRanges(tbl.allDNAMotifs)
      tbl.all.ov <- as.data.frame(findOverlaps(gr.snp, gr.allDNAMotifs, maxgap=shoulder))
      if(nrow(tbl.all.ov) > 0){
         colnames(tbl.all.ov) <- c("snp", "allMotif")
         tbl.allMotifHits <- cbind(tbl.snp[tbl.all.ov$snp,],
                                   tbl.allDNAMotifs[tbl.all.ov$allMotif, c("geneSymbol"), drop=FALSE])
         tbl.allMotifHits <- unique(subset(tbl.allMotifHits, geneSymbol %in% gene.models[[modelName]]$gene))
         if(nrow(tbl.allMotifHits) > 0){
            tbl.allMotifHits$modelName <- modelName
            tbl.allMotifHits$shoulder <- shoulder
            tbl.allMotifHits$source <- "all DNA motifs"
            } # if hits in gene model
         } # if overlap
      } # if allDNA

   tbl.out <- rbind(tbl.fpHits, tbl.dhsMotifHits)
   tbl.out <- rbind(tbl.out, tbl.allMotifHits)
   genes.in.model <- gene.models[[modelName]]$gene
   motif.tfs <- tbl.out$geneSymbol
   match.order <- match(tbl.out$geneSymbol, gene.models[[modelName]]$gene)
   tbl.out$tfRank <- match.order
   tbl.out$pearson <- gene.models[[modelName]]$pearson.coef[match.order]
   rownames(tbl.out) <- NULL

   tbl.out

} # findVariantsInModel
#------------------------------------------------------------------------------------------------------------------------
assay.single.snp <- function(tbl.snp.tf, r, target.gene, tf, tbl.geno, tbl.pheno, mtx)
{
   browser()
   rsid <- tbl.snp.tf$rsid[r]
   printf("%s %s:%s", target.gene, rsid, tf)
   snp.chrom <- tbl.snp.tf$chrom[r]
   snp.loc   <- tbl.snp.tf$pos[r]
   tbl.geno.sub <- subset(tbl.geno, chrom==snp.chrom & start==snp.loc)

   mtx.samples <- colnames(mtx)
   wt.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==0)]
   wt.samples <- intersect(wt.samples, mtx.samples)
   het.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==1)]
   het.samples <- intersect(het.samples, mtx.samples)
   hom.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==2)]
   hom.samples <- intersect(hom.samples, mtx.samples)
   mut.samples <- c(het.samples, hom.samples)

   ad.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis=="AD")$ID)                       # 79
   ctl.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis=="Control")$ID)                 # 73
   pathAging.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis=="Pathologic Aging")$ID)  # 30
   psp.samples <- sub("_TCX", "", subset(tbl.pheno, Diagnosis=="PSP")$ID)                     # 81

   ad.samples <- intersect(ad.samples, mtx.samples)                                           # 79
   ctl.samples <- intersect(ctl.samples, mtx.samples)                                         # 73
   pathAging.samples <- intersect(pathAging.samples, mtx.samples)                             # 30
   psp.samples <- intersect(psp.samples, mtx.samples)                                         # 81

   tbl.plot <- data.frame(target=mtx[target.gene,], tf=mtx[tf,], geno="black", pheno=1, stringsAsFactors=FALSE)
   rownames(tbl.plot) <- colnames(mtx)

   tbl.plot[ctl.samples, "pheno"] <- 15
   tbl.plot[ad.samples, "pheno"]  <- 16
   tbl.plot[pathAging.samples, "pheno"]  <- 17
   tbl.plot[psp.samples, "pheno"]  <- 18

   tbl.plot[wt.samples, "geno"] <- "darkgreen"
   tbl.plot[het.samples, "geno"] <- "pink"
   tbl.plot[hom.samples, "geno"] <- "red"

   with(tbl.plot, plot(tf, target, col=geno, pch=pheno))

   cor(mtx[tf, wt.samples],  mtx[target.gene, wt.samples])
   cor(mtx[tf, het.samples], mtx[target.gene, het.samples])
   cor(mtx[tf, hom.samples], mtx[target.gene, hom.samples])
   cor(mtx[tf, mut.samples], mtx[target.gene, mut.samples])

   cor(mtx[tf, ctl.samples],  mtx[target.gene, ctl.samples])
   cor(mtx[tf, ad.samples], mtx[target.gene, ad.samples])
   cor(mtx[tf, pathAging.samples], mtx[target.gene, pathAging.samples])
   cor(mtx[tf, psp.samples], mtx[target.gene, psp.samples])

   ctl.wt.samples <- intersect(ctl.samples, wt.samples)
   ctl.mut.samples <- intersect(ctl.samples, mut.samples)
   ctl.het.samples <- intersect(ctl.samples, het.samples)
   ctl.hom.samples <- intersect(ctl.samples, hom.samples)

   length(ctl.wt.samples)      # 12
   length(ctl.mut.samples)     # 61
   length(ctl.het.samples)     # 33
   length(ctl.hom.samples)     # 28

   ad.wt.samples <- intersect(ad.samples, wt.samples)
   ad.mut.samples <- intersect(ad.samples, mut.samples)
   ad.het.samples <- intersect(ad.samples, het.samples)
   ad.hom.samples <- intersect(ad.samples, hom.samples)

   length(ad.wt.samples)         # 12
   length(ad.mut.samples)        # 67
   length(ad.het.samples)        # 36
   length(ad.hom.samples)        # 31

   pathAging.wt.samples <- intersect(pathAging.samples, wt.samples)      #  4
   pathAging.mut.samples <- intersect(pathAging.samples, mut.samples)    # 26
   pathAging.het.samples <- intersect(pathAging.samples, het.samples)    # 21
   pathAging.hom.samples <- intersect(pathAging.samples, hom.samples)    #  5

   length(pathAging.wt.samples)        #  4
   length(pathAging.mut.samples)       # 26
   length(pathAging.het.samples)       # 21
   length(pathAging.hom.samples)       #  5

   psp.wt.samples <- intersect(psp.samples, wt.samples)
   psp.mut.samples <- intersect(psp.samples, mut.samples)
   psp.het.samples <- intersect(psp.samples, het.samples)
   psp.hom.samples <- intersect(psp.samples, hom.samples)

   length(psp.wt.samples)             # 18
   length(psp.mut.samples)            # 63
   length(psp.het.samples)            # 36
   length(psp.hom.samples)            # 27

   cor(mtx[tf, ctl.wt.samples],  mtx[target.gene, ctl.wt.samples])                    # 0.49
   cor(mtx[tf, ctl.mut.samples], mtx[target.gene, ctl.mut.samples])                   # 0.68
   cor(mtx[tf, ctl.het.samples],  mtx[target.gene, ctl.het.samples])                  # 0.66
   cor(mtx[tf, ctl.hom.samples], mtx[target.gene, ctl.hom.samples])                   # 0.72

   cor(mtx[tf, ad.wt.samples], mtx[target.gene, ad.wt.samples])                       # 0.69
   cor(mtx[tf, ad.mut.samples], mtx[target.gene, ad.mut.samples])                     # 0.68
   cor(mtx[tf, ad.het.samples], mtx[target.gene, ad.het.samples])                     # 0.74
   cor(mtx[tf, ad.hom.samples], mtx[target.gene, ad.hom.samples])                     # 0.62

   cor(mtx[tf, pathAging.wt.samples], mtx[target.gene, pathAging.wt.samples])         # 0.89
   cor(mtx[tf, pathAging.mut.samples], mtx[target.gene, pathAging.mut.samples])       # 0.36
   cor(mtx[tf, pathAging.het.samples], mtx[target.gene, pathAging.het.samples])       # 0.39
   cor(mtx[tf, pathAging.hom.samples], mtx[target.gene, pathAging.hom.samples])       # 0.02

   cor(mtx[tf, psp.wt.samples], mtx[target.gene, psp.wt.samples])                     # 0.56
   cor(mtx[tf, psp.mut.samples], mtx[target.gene, psp.mut.samples])                   # 0.47
   cor(mtx[tf, psp.het.samples], mtx[target.gene, psp.het.samples])                   # 0.54
   cor(mtx[tf, psp.hom.samples], mtx[target.gene, psp.hom.samples])                   # 0.44




   vec.tf <- mtx[tf,]
   vec.target <- mtx[target.gene,]

   printf("genotype for samples in mtx: %d wt %d het %d hom", length(wt.samples), length(het.samples), length(hom.samples))
   vec.tf.wt       <- as.numeric(mtx[tf, wt.samples])
   vec.tf.mut      <- mtx[tf, c(het.samples, hom.samples)]
   vec.target.wt   <- mtx[target.gene, wt.samples]
   vec.target.mut  <- mtx[target.gene, c(het.samples, hom.samples)]
   correlation <- cor(vec.tf.wt, vec.target.wt)
   title <- sprintf("%s-%s vs %s wt: %5.3f", tf, rsid, target.gene, correlation)
   plot(vec.tf.wt,  vec.target.wt, ylim=c(0,8), xlim=c(0,8), main=title)
   correlation <- cor(vec.tf.mut, vec.target.mut)
   title <- sprintf("%s vs %s mut: %5.3f", tf, target.gene, cor(vec.tf.mut, vec.target.mut))
   plot(vec.tf.mut, vec.target.mut, ylim=c(0,8), xlim=c(0,8), main=title)
   boxplot(vec.tf.wt, vec.tf.mut, vec.target.wt, vec.target.mut)

} # assay.single.snp
#------------------------------------------------------------------------------------------------------------------------
test_assay.single.snp <- function()
{
   printf("--- test_assay.single.snp")
   if(!exists("tbl.mef2c.eqtl.foxp1.motif.zeroShoulder"))   # rs13158247 (score 1.57, IGAP_Pvalue 0.0271) rs244761 (0.42, 0.381), should be good contrast
      print(load("../data/tbl.eqtl.mef2c.2snps.FOXP1.motif.RData"))

   assay.single.snp(tbl.mef2c.eqtl.foxp1.motif.zeroShoulder, 1,  "MEF2C", "FOXP1", tbl.geno, tbl.pheno, mtx.tcx.normalized)

} # test_assay.single.snp
#------------------------------------------------------------------------------------------------------------------------
assay <- function(tbl.snp.tf, target.gene, tf, tbl.geno, tbl.pheno, mtx)
{
   browser()
   stopifnot(tf %in% tbl.snp.tf$geneSymbol)
   tbl.snp.locs <- subset(tbl.snp.tf, geneSymbol==tf)
   all.samples <- colnames(mtx)
   par(mfrow=c(1,3))
   for(r in 1:nrow(tbl.snp.locs)){
      rsid <- tbl.snp.locs$rsid[r]
      printf("%s %s %d/%d", tf, rsid, r, nrow(tbl.snp.locs))
      snp.chrom <- tbl.snp.locs$chrom[r]
      snp.loc   <- tbl.snp.locs$pos[r]
      tbl.geno.sub <- subset(tbl.geno, chrom==snp.chrom & start==snp.loc)
      if(nrow(tbl.geno.sub) == 0)
         next;
      wt.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==0)]
      wt.samples <- intersect(wt.samples, all.samples)
      het.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==1)]
      het.samples <- intersect(het.samples, all.samples)
      hom.samples <- colnames(tbl.geno)[which(subset(tbl.geno, chrom==snp.chrom & start==snp.loc)==2)]
      hom.samples <- intersect(hom.samples, all.samples)
      vec.tf.wt       <- as.numeric(mtx[tf, wt.samples])
      vec.tf.mut      <- mtx[tf, c(het.samples, hom.samples)]
      vec.target.wt   <- mtx[target.gene, wt.samples]
      vec.target.mut  <- mtx[target.gene, c(het.samples, hom.samples)]
      correlation <- cor(vec.tf.wt, vec.target.wt)
      title <- sprintf("%s-%s vs %s wt: %5.3f", tf, rsid, target.gene, correlation)
      plot(vec.tf.wt,  vec.target.wt, ylim=c(0,8), xlim=c(0,8), main=title)
      correlation <- cor(vec.tf.mut, vec.target.mut)
      title <- sprintf("%s vs %s mut: %5.3f", tf, target.gene, cor(vec.tf.mut, vec.target.mut))
      plot(vec.tf.mut, vec.target.mut, ylim=c(0,8), xlim=c(0,8), main=title)
      boxplot(vec.tf.wt, vec.tf.mut, vec.target.wt, vec.target.mut)
      Sys.sleep(5)
      }

} # assay
#------------------------------------------------------------------------------------------------------------------------
test_assay <- function(tbl.snp.tf, tf, tbl.geno, mtx)
{
    assay(tbl.20, "EMX1", tbl.geno, tbl.pheno, mtx.tcx.normalized)
    # EMX1 rs661311 11/27 may have a signal
    # > subset(tbl.20, rsid=="rs661311")
    #     rsid chrom      pos    score iupac ref alt A1 CER_Beta   CER_P TX_Beta  TX_P IGAP_A1 IGAP_OR IGAP_Pvalue RegulomeDB Rsquared_rs254776 Dprime_rs254776 geneSymbol modelName shoulder         source tfRank   pearson
    # rs661311  chr5 88739050 1.872895     Y   C   T  T    -0.11 1.9e-06   -0.03 0.348       T    0.95      0.0134          7             0.485               1       EMX1 mef2c.tcx       20 all DNA motifs      5 0.7821005


} # assay
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) {

   context = init.context()
   socket = init.socket(context,"ZMQ_REP")
   bind.socket(socket, sprintf("tcp://*:%d", PORT))

   errorFunction <- function(condition){
     printf("==== exception caught ===")
     print(as.character(condition))
     response <- list(cmd="handleError", status="error", callback="", payload=as.character(condition));
     send.raw.string(socket, toJSON(response))
     };

   while(TRUE) {
      tryCatch({
        printf("top of receive/send loop")
        raw.message <- receive.string(socket)
        msg = fromJSON(raw.message)
        printf("--- msg:")
        print(msg)

        stopifnot(is.list(msg))
        expected.fields <- c("cmd", "status", "callback", "payload")
        stopifnot(all(expected.fields %in% names(msg)))

        response <- processWellStructuredMessage(msg)
        json.string <- toJSON(response, dataframe="values")
        send.raw.string(socket, json.string)
        Sys.sleep(1)
        }, error=errorFunction); # tryCatch
     } # while (TRUE)

} # if !interactive()
#------------------------------------------------------------------------------------------------------------------------
test_getFootprints <- function()
{
   printf("--- test_getFootprints")

    # fp.roi specifies the range for which we already have footprints
   test.roi <- fp.roi
   test.roi$start <- test.roi$start + 10000
   test.roi$end <- test.roi$start + 20000
   test.roi.string <- with(test.roi, sprintf("%s:%d-%d", chrom, start, end))
   fp.roi.string <- with(fp.roi, sprintf("%s:%d-%d", chrom, start, end))

   tbl.test <- getFootprints(test.roi.string)
   checkEquals(nrow(tbl.test), 203)
   checkEquals(ncol(tbl.test), 12)
   checkEquals(colnames(tbl.test)[1:3], c("chrom", "start", "end"))

   outOfRange.roi <- sprintf("%s:%d-%d", fp.roi$chrom, fp.roi$end + 10, fp.roi$end+1000)
   tbl.outOfRange <- getFootprints(outOfRange.roi)
   checkTrue(nrow(tbl.outOfRange) > 20)
   checkEquals(ncol(tbl.outOfRange), 12)
   checkEquals(colnames(tbl.outOfRange)[1:3], c("chrom", "start", "end"))

} # test_getFootprints
#------------------------------------------------------------------------------------------------------------------------
test_getDHSRegions <- function()
{
   printf("--- test_getDHSRegions")

   test.roi <- fp.roi
   test.roi$start <- test.roi$start + 10000
   test.roi$end <- test.roi$start + 20000
   test.roi.string <- with(test.roi, sprintf("%s:%d-%d", chrom, start, end))
   fp.roi.string <- with(fp.roi, sprintf("%s:%d-%d", chrom, start, end))

   tbl.dhs.new <- getDHSRegions(fp.roi.string)
   checkEquals(dim(tbl.dhs), dim(tbl.dhs.new))
   checkEquals(colnames(tbl.dhs)[1:3], c("chrom", "start", "end"))

} # test_getDHSRegions
#------------------------------------------------------------------------------------------------------------------------
test_getDHSMotifs <- function()
{
   printf("--- test_getDHSMotifs")

   test.roi <- fp.roi
   test.roi$start <- test.roi$start + 10000
   test.roi$end <- test.roi$start + 20000
   test.roi.string <- with(test.roi, sprintf("%s:%d-%d", chrom, start, end))
   fp.roi.string <- with(fp.roi, sprintf("%s:%d-%d", chrom, start, end))

   tbl.dhs.new <- getDHSMotifs(fp.roi.string)
   checkEquals(dim(tbl.dhsMotifs), dim(tbl.dhs.new))
   checkEquals(colnames(tbl.dhs.new)[1:3], c("chrom", "start", "end"))

} # test_getDHSRegions
#------------------------------------------------------------------------------------------------------------------------
test_findVariantsInModel <- function()
{
   printf("--- test_findVariantsInModel")
   tbl.00.fp <- findVariantsInModel(modelName="mef2c.tcx", shoulder=0, motifSource=c("footprints"))
   checkEquals(dim(tbl.00.fp), c(3, 24))

   tbl.00.dhs <- findVariantsInModel(modelName="mef2c.tcx", shoulder=0, motifSource=c("dhs"))
   checkEquals(dim(tbl.00.dhs), c(2,24))

   tbl.00.allDNA <- findVariantsInModel(modelName="mef2c.tcx", shoulder=0, motifSource=c("allDNA"))
   checkEquals(dim(tbl.00.allDNA), c(14, 24))

   tbl.15 <- findVariantsInModel(modelName="mef2c.tcx", shoulder=15, motifSource=c("footprints", "dhs"))
   checkEquals(dim(tbl.15), c(15, 24))

   tbl.xtab <- as.data.frame(table(tbl.00.fp$source), stringsAsFactors=FALSE)
   checkEquals(tbl.xtab$Var1, c("hint+wellington, 16+20 footprints"))
   checkEquals(tbl.xtab$Freq, c(3))

   tbl.xtab <- as.data.frame(table(tbl.00.fp$geneSymbol), stringsAsFactors=FALSE)
   checkEquals(tbl.xtab$Var1, c("FOXP1", "MEF2A"))
   checkEquals(tbl.xtab$Freq, c(2, 1))

      # for a digestible view:
   coi <- c(1, 2, 3, 4, 6, 7, 17, 19, 23, 24, 20, 21, 22)
      #    tbl.00[, coi]
   tbl.05 <- findVariantsInModel(modelName="mef2c.tcx", motifSources=c("footprints", "dhs", "allDNA"),  shoulder=5)
   checkEquals(dim(tbl.05), c(34, 24))
   tbl.xtab <- as.data.frame(table(tbl.05$source), stringsAsFactors=FALSE)
   checkEquals(tbl.xtab$Var1, c("all DNA motifs", "encode DHS motifs", "hint+wellington, 16+20 footprints"))
   checkEquals(tbl.xtab$Freq, c(25, 3, 6))

   tbl.xtab <- as.data.frame(table(tbl.05$geneSymbol), stringsAsFactors=FALSE)
   checkEquals(tbl.xtab$Var1, c("BHLHE22","EMX1", "FOXP1", "FOXP2", "HLF", "LBX2", "MEF2A", "NFE2L2", "SP3"))
   checkEquals(tbl.xtab$Freq, c(1, 10, 2, 3, 2, 10, 3, 2, 1))

} # test_findVariantsInModel
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getFootprints()
   test_getDHSRegions()
   test_getDHSMotifs()
   test_findVariantsInModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
