library(rzmq)
library(jsonlite)
library(trena)
PORT <- 5547
#------------------------------------------------------------------------------------------------------------------------
# load expression matrices on startup
#
matrix.files <- c("data/mtx.protectedAndExposed.RData",
                  "data/gtex.fibroblast.RData",
                  "data/gtex.primary.RData")

expression.matrix.files <- list(protectedAndExposed="data/mtx.protectedAndExposed.RData",
                                gtexFibroblast="data/gtex.fibroblast.RData",
                                gtexPrimary="data/gtex.primary.RData")
expression.matrices <- list()

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
   else if(msg$cmd == "getFootprintsInRegion"){
      stopifnot("roi" %in% names(msg$payload))
      roi <- msg$payload$roi
      tbl.fp.as.list <- getFootprints(roi)
      response <- list(cmd=msg$callback, status="success", callback="", payload=tbl.fp.as.list);
      }
   else if(msg$cmd == "listSharedData"){
      printf("--- about to get shared data filenames")
      print(list.files("/home/trena/work"))
      print(list.files("/home/trena/work/shared"))
      filenames <- list.files("/home/trena/work/shared")
      response <- list(cmd=msg$callback, status="success", callback="", payload=filenames);
      }
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
   print(6)
   dataFrameToPandasFriendlyList(tbl.reg)

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
