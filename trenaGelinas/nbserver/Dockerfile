FROM ipytrenaviz:latest
USER root
WORKDIR /home/jovyan
ADD . /home/jovyan/trenaViz/
RUN chown -R jovyan:users /home/jovyan
USER jovyan
WORKDIR /home/jovyan/trenaViz
ADD custom.css /home/jovyan/.jupyter/custom/

