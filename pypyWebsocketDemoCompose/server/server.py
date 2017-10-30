import zmq
import time
import sys

port = "5557"

context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:%s" % port)

while True:
    print("top of server recv loop");
    message = socket.recv_string()
    print("Received request: %s" % message)
    time.sleep (1)
    socket.send_string(message.upper())
