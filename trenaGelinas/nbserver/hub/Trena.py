import zmq, json
import pandas as pd
from ipyTrenaViz import *
import time

class Trena:

    def __init__(self, genomeName):
       socketContext = zmq.Context();
       self.trenaServer = socketContext.socket(zmq.REQ)
       self.trenaServer.connect("tcp://trena:%s" % "5547")
       self.tv = ipyTrenaViz()
       display(self.tv)
       self.tv.setGenome(genomeName)

    def display(self):
       display(self.tv)

    def ping(self):
        msg = {'cmd': 'ping', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        return(response)

    def showGenomicRegion(self, regionString):
        self.tv.showGenomicRegion(regionString);

    def getGenomicRegion(self):
        return(self.tv.getBrowserState()["chromLocString"]);

    def dataFrameFrom3partList(self, list):
        data = list['tbl']
        rownames = list['rownames']
        colnames = list['colnames']
        df = pd.DataFrame(data)
        df.columns = colnames
        rownameList = {}
        for i in range(len(rownames)):
          rownameList[i] = rownames[i]
        df = df.rename(rownameList)
        return(df)

    def getExpressionMatrixNames(self):
        msg = {'cmd': 'getExpressionMatrixNames', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload)

    def summarizeExpressionMatrices(self):
        msg = {'cmd': 'summarizeExpressionMatrices', 'status': 'request', 'callback': '', 'payload': ''}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(self.dataFrameFrom3partList(payload))

    def getFootprintsInRegion(self, display):
        payload = {"roi": self.getGenomicRegion()}
        msg = {'cmd': 'getFootprintsInRegion', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        tbl = self.dataFrameFrom3partList(payload)
        if(display):
           self.tv.addBedTrackFromDataFrame(tbl)
        return(tbl)

    def displayFootprints(self, url):
        self.tv.addBedTrackFromDataFrame(url)

    def addBedTrackFromDataFrame(self, tbl, trackName, trackMode, color, trackHeight=200):
        return(self.tv.addBedTrackFromDataFrame(tbl, trackName, trackMode, color, trackHeight))

    def displayGraph(self, filename, modelName):
        self.tv.displayGraph(filename, modelName)

    def sessionInfo(self):
        msg = {'cmd': 'getSessionInfo', 'status': 'request', 'callback': '', 'payload': ""}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload);

    def listSharedData(self):
        msg = {'cmd': 'listSharedData', 'status': 'request', 'callback': '', 'payload': ""}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(payload);

    def createGeneModel(self, targetGene,  matrixName):
        payload = {'targetGene': targetGene,
                   'matrixName': matrixName}
        msg = {'cmd': 'createGeneModel', 'status': 'request', 'callback': '', 'payload': payload}
        self.trenaServer.send_string(json.dumps(msg))
        response = json.loads(self.trenaServer.recv_string())
        payload = response["payload"]
        return(self.dataFrameFrom3partList(payload))


