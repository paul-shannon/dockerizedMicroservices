FROM pshannon/trenaprereqs:latest
user root
WORKDIR /home/trena
COPY ./installTrena.sh .
COPY ./.Rprofile .
COPY ./ .
user trena
RUN bash ./installTrena.sh


