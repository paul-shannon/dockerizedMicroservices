FROM pshannon/trenaserver:latest
USER root
WORKDIR /home/trena
RUN mkdir privateData
COPY myBashCustomizations /home/trena
RUN cat myBashCustomizations >> /home/trena/.bashrc
COPY server.R /home/trena/work/
COPY trena_1.1.5.tar.gz /home/trena/work/
COPY data/* /home/trena/privateData/
USER root
RUN mkdir sharedData
RUN chmod 777 sharedData
RUN chown -R trena:users /home/trena/work
RUN chown -R trena:users /home/trena/privateData
RUN chown -R trena:users /home/trena/sharedData
USER trena
RUN R CMD INSTALL work/trena_1.1.5.tar.gz
WORKDIR /home/trena/work
CMD ["R", "-f", "server.R"]
