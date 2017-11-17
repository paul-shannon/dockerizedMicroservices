FROM trenaserver:latest
WORKDIR /home/trena/work
ADD server.R /home/trena/work/
COPY trena_1.1.5.tar.gz /home/trena/work/
ADD data/ /home/trena/work/data/
USER root
RUN chown -R trena:users /home/trena/work
USER trena
RUN R CMD INSTALL trena_1.1.5.tar.gz
RUN mkdir shared
RUN chmod 777 shared
CMD ["R", "-f", "server.R"]
