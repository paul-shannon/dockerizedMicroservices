FROM bioconductor/devel_core2:latest


COPY . /usr/local/src/myscripts2
WORKDIR /usr/local/src/myscripts2
RUN R CMD INSTALL /usr/local/src/myscripts2/Demo_0.99.0.tar.gz

#RUN echo "biocGet <- function(pkgName){source("http://bioconductor.org/biocLite.R"); biocLite(pkgName)}" > ~/.Rprofile
##RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" >> ~/.Rprofile
#RUN Rscript -e "install.packages('glmnet')"
##RUN Rscript -e "install.packages('glmnet')"
##RUN Rscript -e "install.packages('lassopv')"
##RUN Rscript -e "install.packages('randomForest')"
##RUN Rscript -e "install.packages('flare')"
##RUN Rscript -e "install.packages('vbsr')"
##RUN Rscript -e "install.packages('RPostgreSQL')"
###RUN Rscript -e "install.packages('RMySQL')"
##RUN Rscript -e "install.packages('BSgenome.Hsapiens.UCSC.hg38')"
##RUN Rscript -e "install.packages('BSgenome.Hsapiens.UCSC.hg19')"
##RUN Rscript -e "install.packages('BSgenome.Mmusculus.UCSC.mm10')"
##RUN Rscript -e "install.packages('org.Hs.eg.db')"
##RUN Rscript -e "install.packages('SNPlocs.Hsapiens.dbSNP144.GRCh38')"
##RUN R CMD INSTALL BrowserViz_1.9.15.tar.gz
##RUN R CMD INSTALL Demo_0.99.0.tar.gz
##RUN R CMD INSTALL splitstackshape_1.4.2.tar.gz
##RUN R CMD INSTALL MotifDb_1.19.17.tar.gz
###RUN R CMD INSTALL trena_0.99.187.tar.gz
###RUN R CMD INSTALL trenaViz_0.99.21.tar.gz
###RUN R CMD INSTALL trenaGWAS01_0.99.6.tar.gz
#
#USER trenaUser
#WORKDIR /home/trenaUser

CMD ["Rscript", "myscript.R"]
