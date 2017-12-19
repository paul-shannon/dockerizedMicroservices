FROM pshannon/trenabase:latest
USER root
RUN apt-get install libzmq3-dev -y
RUN apt-get install libpq-dev -y
USER trena
WORKDIR /home/trena
COPY . /home/trena/
RUN R -f installTrenaServerPackages.R
RUN R CMD INSTALL MotifDb_1.19.18.tar.gz
RUN R CMD INSTALL trena_1.1.10.tar.gz


