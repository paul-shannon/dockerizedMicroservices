library(trenaViz)
PORT.RANGE <- 8000:8020
if(!exists("tv")){
   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
   setGenome(tv, "hg38")
   }

tbl.fpLeft <- read.table("fpLeft.tsv", sep="\t", as.is=TRUE, header=TRUE)
