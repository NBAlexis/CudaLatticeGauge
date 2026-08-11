import numpy as np
from matplotlib import pyplot as plt

from AutoCorrelation import AutoCorrelationSingleVariable
from Distance.UsefulFunctions import findDistances, findK
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMatlabArray

folderstart = "G:\\"

klst = findK(folderstart, "distancer")
distlst = findDistances(folderstart, "distancer")

betalst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]
col1 = [3*i for i in range(len(distlst))]
col2 = [3*i+1 for i in range(len(distlst))]
col3 = [3*i+2 for i in range(len(distlst))]

k1lst = []
k2lst = []
k3lst = []
k4lst = []
k5lst = []
k6lst = []

for i in range(len(betalst)):
    wilsonloop = LoadMathematicaCSV(f"{folderstart}QuenchWilsonLoop\\DistanceR\\BSQ_DistanceR__{betalst[i]}_wilsonloops.csv")
    # every 3 cols correspond to one result
    wilsonloop2 = np.real(wilsonloop[:, col1] + wilsonloop[:, col2] + wilsonloop[:, col3])/3
    for k in range(6):
        rk = distlst[klst == k + 1]
        datak = wilsonloop2[:, klst == k + 1]
        # print(np.shape(datak))
        argidx = rk.argsort()
        rk = rk[argidx]
        datak = datak[:, argidx]
        rlst = []
        datalst = []
        lastr = 0
        for j in range(1, len(rk)):
            # check whether j is as same as j - 1
            # if not, gather from lastr to j - 1
            if rk[j] != rk[j - 1]:
                rlst.append(rk[lastr])
                if lastr == j - 1:
                    datalst.append(datak[:, lastr])
                else:
                    # print(f'k={k}, lastr={lastr}, j={j}, rk={rk[j]}')
                    # print(np.shape(datak[:, lastr:j]))
                    # print(np.shape(np.sum(datak[:, lastr:j], axis=1)))
                    datalst.append(np.sum(datak[:, lastr:j], axis=1) / (j - lastr))
                lastr = j
        rlst.append(rk[lastr])
        if lastr != len(rk) - 1:
            datalst.append(np.sum(datak[:, lastr:len(rk)], axis=1) / (len(rk) - lastr))
        else:
            datalst.append(datak[:, lastr])
        rlst = np.array(rlst)
        datalst = np.array(datalst)
        rlst = np.sqrt(rlst)
        vs = []
        ss = []
        ts = []
        for j in range(len(datalst)):
            v, s = JacknifeMean(datalst[j], f'b:{i}-beta:{betalst[i]}-r:{k}-{j}')
            _, _, t = AutoCorrelationSingleVariable(datalst[j])
            vs.append(v)
            ss.append(s * np.sqrt(2 * t))
            ts.append(t)
        vs = np.array(vs)
        ss = np.array(ss)
        ts = np.array(ts)
        if 0 == k:
            k1lst.append(rlst)
            k1lst.append(vs)
            k1lst.append(ss)
            k1lst.append(ts)
        elif 1 == k:
            k2lst.append(rlst)
            k2lst.append(vs)
            k2lst.append(ss)
            k2lst.append(ts)
        elif 2 == k:
            k3lst.append(rlst)
            k3lst.append(vs)
            k3lst.append(ss)
            k3lst.append(ts)
        elif 3 == k:
            k4lst.append(rlst)
            k4lst.append(vs)
            k4lst.append(ss)
            k4lst.append(ts)
        elif 4 == k:
            k5lst.append(rlst)
            k5lst.append(vs)
            k5lst.append(ss)
            k5lst.append(ts)
        elif 5 == k:
            k6lst.append(rlst)
            k6lst.append(vs)
            k6lst.append(ss)
            k6lst.append(ts)

k1lst = np.array(k1lst)
k2lst = np.array(k2lst)
k3lst = np.array(k3lst)
k4lst = np.array(k4lst)
k5lst = np.array(k5lst)
k6lst = np.array(k6lst)

np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK1.csv', k1lst, delimiter=',')
np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK2.csv', k2lst, delimiter=',')
np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK3.csv', k3lst, delimiter=',')
np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK4.csv', k4lst, delimiter=',')
np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK5.csv', k5lst, delimiter=',')
np.savetxt(f'{folderstart}\\WilsonLoopRes\\BSQDistanceRK6.csv', k6lst, delimiter=',')
