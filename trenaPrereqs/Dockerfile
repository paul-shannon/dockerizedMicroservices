FROM bioconductor/devel_core2:latest
RUN useradd -ms /bin/bash trena
USER root
RUN apt-get -y install libmariadb-client-lgpl-dev
RUN apt-get install libpq-dev -y
user trena
WORKDIR /home/trena
RUN mkdir /home/trena/work
RUN mkdir /home/trena/library
COPY . /home/trena/
RUN cat /home/trena/myBashCustomizations >> /home/trena/.bashrc
RUN R CMD INSTALL -l ~/library SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz
RUN rm SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz
RUN R -f installTrenaPrequisites.R
