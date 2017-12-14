library(rzmq)
library(jsonlite)
library(trena)
library(trenaViz) # used here only for buildMultiModelGraph and addGeneLayout
PORT <- 5547
#------------------------------------------------------------------------------------------------------------------------
# load expression matrices on startup
#
expression.matrix.files <- list(protectedAndExposed="../privateData/mtx.protectedAndExposedLowVarianceGenesEliminated.RData",
                                gtexFibroblast="../privateData/gtex.fibroblast.RData",
                                gtexPrimary="../privateData/gtex.primary.RData")
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
getFootprints <- function(roi)
{
   trena <- Trena("hg38")
   source.1 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_16"
   source.2 <- "postgres://bddsrds.globusgenomics.org/skin_wellington_20"
   source.3 <- "postgres://bddsrds.globusgenomics.org/skin_hint_16"
   source.4 <- "postgres://bddsrds.globusgenomics.org/skin_hint_20"
   sources <- c(source.1, source.2, source.3, source.4)
   names(sources) <- c("well_16", "well_20", "hint_16", "hint_20")

   sources <- sources[4];

   targetGene <- "COL1A1"
   tss <- 50201632

   cls <- parseChromLocString(roi)   # a trena function

   x <- getRegulatoryChromosomalRegions(trena, cls$chrom, cls$start, cls$end, sources, targetGene, tss)
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

   printf("---- getFootprints, after associateTranscriptionFactors")
   print(head(tbl.reg))

   tbl.bed <- tbl.reg[, c("chrom", "motifStart", "motifEnd")]
   tbl.bed$fpName <- paste(tbl.reg$motifName, tbl.reg$db, sep="_")
      # nasty hack, working around igv.js, npm, ipywidgets, docker, ???: last line in tbl is dropped
   tbl.bed <- rbind(tbl.bed, tbl.bed[nrow(tbl.bed),])
   #printf("--- writing tbl.bed")
   #write.table(tbl.bed, file="/home/trena/sharedData/tbl.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
   #printf("--- after write")
   print(6)
   tbl.reg

} # getFootprints
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
