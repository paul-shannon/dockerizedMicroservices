import sys
sys.path.append("../hub")
from trena import *

def runTests():
    testConstructor()


def testConstructor():
    print("--- testConstuctor")
    name = "ABCD"
    trena = Trena(name)
    assert(trena.getName() == name);

if(__name__ == "__main__"):
    runTests()


