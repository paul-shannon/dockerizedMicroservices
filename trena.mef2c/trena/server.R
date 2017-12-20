library(rzmq)
library(jsonlite)
library(trena)
library(trenaViz) # used here only for buildMultiModelGraph and addGeneLayout
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
for(model.name in mef2c.model.names){
   gene.models[[model.name]] <- eval(parse(text=model.name))
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
#     dim(tbl.dhsMotifs) # 28872 9
#     save(tbl.dhsMotifs, file="../data/tbl.dhsMotifs.RData")
#------------------------------------------------------------------------------------------------------------------------
# load the open chromatin (dhs) table
printf("loading dhs regions: %s", load("../data/tbl.dhs.RData"))
printf("loading motifs in dhs regions: %s", load("../data/tbl.dhsMotifs.RData"))
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
      tbl.dhs.roi <- getDHSRegions(roi)
      if(nrow(tbl.dhs.ro1) == 0){
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
      print(head(tbl.motifs.tfs))

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
runTests <- function()
{
   test_getFootprints()
   test_getDHSRegions()
   test_getDHSMotifs()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
