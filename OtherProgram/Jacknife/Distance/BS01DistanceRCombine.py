import numpy as np
from matplotlib import pyplot as plt

from AutoCorrelation import AutoCorrelationSingleVariable
from Distance.UsefulFunctions import findDistances, findK
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMatlabArray

folderstart = "G:\\"

klst = findK(folderstart, "distancer")
distlst = findDistances(folderstart, "distancer")

betalst = ["530", "535", "540", "545", "550", "555", "560", "565", "570", "575", "580"]
col1 = [3*i for i in range(len(distlst))]
col2 = [3*i+1 for i in range(len(distlst))]
col3 = [3*i+2 for i in range(len(distlst))]

alllst = []

for i in range(len(betalst)):
    wilsonloop = LoadMathematicaCSV(f"{folderstart}BS01Wilsonloops\\DistanceR\\BS01_DistanceR__{betalst[i]}_wilsonloops.csv")
    # every 3 cols correspond to one result
    wilsonloop2 = np.real(wilsonloop[:, col1] + wilsonloop[:, col2] + wilsonloop[:, col3])/3
    vs = []
    ss = []
    ts = []
    for k in range(6):
        rk = distlst[klst == k + 1]
        datak = wilsonloop2[:, klst == k + 1]
        # print(np.shape(datak))
        argidx = rk.argsort()
        rk = rk[argidx]
        datak = datak[:, argidx]
        count = 0
        datalst = None
        for j in range(0, len(rk)):
            # check whether j is as same as j - 1
            # if not, gather from lastr to j - 1
            if rk[j] >= 36:
                if datalst is None:
                    datalst = datak[:, j]
                    count = count + 1
                else:
                    datalst = datalst + datak[:, j]
                    count = count + 1
        datalst = datalst/count
        v, s = JacknifeMean(datalst, f'b:{i}-beta:{betalst[i]}-r:{k}')
        _, _, t = AutoCorrelationSingleVariable(datalst)
        vs.append(v)
        ss.append(s * np.sqrt(2 * t))
        ts.append(t)
    alllst.append(vs)
    alllst.append(ss)
    alllst.append(ts)

np.savetxt(f'{folderstart}\\WilsonLoopRes\\BS01DistanceRCombine.csv', alllst, delimiter=',')
