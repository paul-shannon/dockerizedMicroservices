import zmq, json
from ipyTrenaViz import *

class Trena:

    def __init__(self, name):
       self.name = name;
       socketContext = zmq.Context();
       self.trenaServer = socketContext.socket(zmq.REQ)
       self.trenaServer.connect("tcp://trena:%s" % "5547")
       self.tv = ipyTrenaViz()

    def display(self):
       display(self.tv)
       self.tv.setGenome()

    def getName(self):
        return self.name

    def ping(self):
        msg = {'cmd': 'ping', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        return(response)




