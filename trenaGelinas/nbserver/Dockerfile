FROM pshannon/ipytrenaviz:latest
USER root
WORKDIR /home/jovyan
COPY myBashCustomizations /home/jovyan/
RUN cat myBashCustomizations >> /home/jovyan/.bashrc

COPY col1a1.ipynb /home/jovyan/trenaViz/
RUN mkdir /home/jovyan/trenaViz/hub
COPY hub/* /home/jovyan/trenaViz/hub/

# ADD . /home/jovyan/trenaViz/
RUN mkdir /home/jovyan/trenaViz/privateData
RUN mkdir /home/jovyan/trenaViz/sharedData
COPY g-5nodes4edges.json /home/jovyan/trenaViz/sharedData
COPY tmp.bed /home/jovyan/trenaViz/sharedData

RUN chown -R jovyan:users /home/jovyan

USER jovyan
WORKDIR /home/jovyan/trenaViz
ADD custom.css /home/jovyan/.jupyter/custom/
