library(rzmq)
library(jsonlite)
context = init.context()
socket = init.socket(context,"ZMQ_REP")
bind.socket(socket,"tcp://*:5557")
while(TRUE) {
   request = receive.string(socket)
   printf("server.R received REQUEST: %s-%s", request, "from server.R")
   response <- sprintf("%s (from pyRWebSocketDemoCompose/server.R)", toupper(request))
   send.raw.string(socket, response)
   Sys.sleep(1)
   }
