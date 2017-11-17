FROM ipytrenaviz:latest
USER root
WORKDIR /home/jovyan
ADD . /home/jovyan/trenaViz/
RUN mkdir /home/jovyan/trenaViz/nbData
RUN chown -R jovyan:users /home/jovyan

USER jovyan
WORKDIR /home/jovyan/trenaViz
ADD custom.css /home/jovyan/.jupyter/custom/
RUN mkdir shared
RUN chmod 777 shared
