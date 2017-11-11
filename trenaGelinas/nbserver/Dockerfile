FROM ipytrenaviz:latest
WORKDIR /home/jovyan
ADD . /home/jovyan/trenaViz/
USER root
RUN chown -R jovyan:users /home/jovyan
USER jovyan
WORKDIR /home/jovyan/trenaViz
ADD custom.css /home/jovyan/.jupyter/custom/

