FROM pshannon/trenaserver:latest
USER root
WORKDIR /home/trena
RUN mkdir data
COPY myBashCustomizations /home/trena
RUN cat myBashCustomizations >> /home/trena/.bashrc
COPY server.R /home/trena/work/
COPY trena_1.1.10.tar.gz /home/trena/work/
COPY trenaViz_0.99.34.tar.gz /home/trena/work
COPY BrowserViz_1.9.15.tar.gz /home/trena/work
COPY data/*.RData /home/trena/data/
COPY data/README.txt /home/trena/data/
USER root
RUN mkdir sharedData
RUN chmod 777 sharedData
RUN chown -R trena:users /home/trena/work
RUN chown -R trena:users /home/trena/data
RUN chown -R trena:users /home/trena/sharedData
USER trena
RUN R CMD INSTALL work/trena_1.1.10.tar.gz
RUN R CMD INSTALL work/BrowserViz_1.9.15.tar.gz
RUN R CMD INSTALL work/trenaViz_0.99.34.tar.gz
WORKDIR /home/trena/work
CMD ["R", "-f", "server.R"]
