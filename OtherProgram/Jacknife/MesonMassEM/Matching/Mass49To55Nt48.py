from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

disk = "G:\\BS01\\Meson\\"
beta = ["49", "50", "51", "52", "53", "54", "55"]
T = 48
d = 2
maxConf = 800

for b in beta:
    allpv = []
    print("beta = {}", b)
    for pv in range(4):
        filename = "{}BS01__mesonsimple{}_{}.csv". format(disk, pv, b)
        testarray = LoadMathematicaCSV(filename)
        testarray = testarray[100:maxConf, :]
        res = JacknifeMeasonSimple(testarray, T, d)
        print("{} {} = {}".format(b, pv, res))
    allpv.append(res)


