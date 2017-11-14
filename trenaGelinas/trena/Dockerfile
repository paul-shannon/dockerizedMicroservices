FROM trenaserver:latest
WORKDIR /home/trena/work
ADD server.R /home/trena/work
ADD data/ /home/trena/work/data/
USER root
RUN chown -R trena:users /home/trena/work
USER trena
CMD ["R", "-f", "server.R"]
