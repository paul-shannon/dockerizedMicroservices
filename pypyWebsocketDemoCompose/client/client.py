import zmq
f = open ("clientOutput.text", "w")
f.write("about to get zmq context\n");
f.flush()
socketContext = zmq.Context()
socket = socketContext.socket(zmq.REQ)
socket.connect("tcp://websocketserver:%s" % '5557')
for i in range(5):
   print("About to ***SEnD*** string to server")
   socket.send_string("hello from python: %d" % i)
   print("RETuRNed from server: %s", socket.recv_string())
