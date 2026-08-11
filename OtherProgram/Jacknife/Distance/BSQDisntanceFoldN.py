import numpy as np
from matplotlib import pyplot as plt

from AutoCorrelation import AutoCorrelationSingleVariable
from Distance.UsefulFunctions import findDistances, findK
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean

folderstart = "G:\\"

klst = findK(folderstart, "distancefoldn")
distlst = findDistances(folderstart, "distancefoldn")
print(klst)
print(distlst)
print(len(distlst))

betalst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]
col1 = [6*i for i in range(len(distlst))]
col2 = [6*i+1 for i in range(len(distlst))]
col3 = [6*i+2 for i in range(len(distlst))]
col4 = [6*i+3 for i in range(len(distlst))]
col5 = [6*i+4 for i in range(len(distlst))]
col6 = [6*i+5 for i in range(len(distlst))]

kalllst = []

for i in range(len(betalst)):
    wilsonloop = LoadMathematicaCSV(f"{folderstart}QuenchWilsonLoop\\DistanceFoldN\\BSQ_DistanceFoldN__{betalst[i]}_wilsonloops.csv")
    # every 3 cols correspond to one result
    wilsonloop2 = np.real(wilsonloop[:, col1] + wilsonloop[:, col2] + wilsonloop[:, col3] + wilsonloop[:, col4] + wilsonloop[:, col5] + wilsonloop[:, col6])/6
    for k in range(6):
        distk = distlst[klst == k + 1]
        loopk = wilsonloop2[:, klst == k + 1]
        vs = []
        ss = []
        ts = []
        for j in range(len(distk)):
            v, s = JacknifeMean(loopk[:, j], f'b:{i}-beta:{betalst[i]}-k:{k}-r:{j}')
            _, _, t = AutoCorrelationSingleVariable(loopk[:, j])
            vs.append(v)
            ss.append(s * np.sqrt(2 * t))
            ts.append(t)
        vs = np.array(vs)
        ss = np.array(ss)
        ts = np.array(ts)
        kalllst.append(distk)
        kalllst.append(vs)
        kalllst.append(ss)
        kalllst.append(ts)

np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceFoldN.csv', kalllst, delimiter=',')

