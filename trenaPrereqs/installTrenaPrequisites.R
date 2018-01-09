printf <- function(...) print(noquote(sprintf(...)))
my.user.library <- "~/library"
if(!my.user.library %in% .libPaths())
    .libPaths(my.user.library)

source("http://bioconductor.org/biocLite.R")

snp.package <- "SNPlocs.Hsapiens.dbSNP150.GRCh38"
snp.package.tarball <- "SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz"


data.pkgs <- c(snp.package,
               "BSgenome",
               "BSgenome.Hsapiens.UCSC.hg38",
               "BSgenome.Hsapiens.UCSC.hg19",
               "BSgenome.Mmusculus.UCSC.mm10",
               "org.Hs.eg.db",
               "org.Mm.eg.db",
               "BSgenome.Athaliana.TAIR.TAIR9"
               )

for(data.pkg in data.pkgs){
   suppressWarnings(
      needed <- !require(data.pkg, character.only=TRUE, lib.loc=my.user.library, quiet=TRUE)
      )
   printf("%s needed? %s", data.pkg, needed)
   if(needed){
     if(grepl("SNPlocs.Hsapiens", data.pkg)){   # docker build has trouble with download
        printf("installing snplocs from local tarball: %s", snp.package.tarball)
        install.packages(snp.package.tarball, lib="~/library", repos=NULL)
        }
     else{
        biocLite(data.pkg,quiet=FALSE, lib=my.user.library)
        }
      }
   } # for



code.pkgs <- c("GenomicRanges",
               "Biostrings",
               "BiocParallel",
               "DBI",
               "RPostgreSQL",
               "RMySQL",
               "RSQLite",
               "glmnet",
               "lassopv",
               "randomForest",
               "flare",
               "vbsr",
               "stringr",
               "httpuv",
               "colorspace",
               "annotate",
               "MotifDb",
               "splitstackshape",
               "RUnit")

for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE, lib.loc=my.user.library, quiet=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocLite(code.pkg, quiet=FALSE, lib=my.user.library)
   } # for




