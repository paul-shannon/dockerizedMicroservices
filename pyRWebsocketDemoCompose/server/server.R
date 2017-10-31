library(rzmq)
library(jsonlite)
context = init.context()
socket = init.socket(context,"ZMQ_REP")
bind.socket(socket,"tcp://*:5557")
while(TRUE) {
   request = receive.string(socket)
   printf("received request: %s", request)
   response <- toupper(request)
   send.raw.string(socket, response)
   Sys.sleep(1)
   }
