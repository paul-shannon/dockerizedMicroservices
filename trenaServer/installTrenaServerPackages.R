printf <- function(...) print(noquote(sprintf(...)))

source("http://bioconductor.org/biocLite.R")

code.pkgs <- c("jsonlite",
               "rzmq",
               "RPostgreSQL",
               "glmnet",
               "RUnit")


biocLite(code.pkgs, lib="/home/trena/library")

data.pkgs <- c()

for(data.pkg in data.pkgs){
   suppressWarnings(
      needed <- !require(data.pkg, character.only=TRUE, quiet=TRUE)
      )
   printf("%s needed? %s", data.pkg, needed)
   if(needed)
      biocLite(data.pkg,quiet=FALSE, lib="/home/trena/library")
   } # for


