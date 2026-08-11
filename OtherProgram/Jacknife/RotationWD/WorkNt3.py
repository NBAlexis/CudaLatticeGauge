import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import JacknifeMean, JacknifeCumulant, LoadMathematicaCSV
from Visualization import errorbar, LineStyle, MarkerStyle

# filenames = "F:\\{}\\Polyakov\\NC{}__polyakov_Nt3_{}.csv"
filenames = "G:\\RotationWD\\Nt3Polyakov\\{}\\NC{}__polyakov_Nt3_{}.csv"
writefolder = "G:\\RotationWD\\ExportCSV\\Nt3\\"
useRealInsteadOfAbs = True
folders = ["46916", "47216", "47516", "47816", "48116", "48416", "48716", "49016", "49316", "49616", "49916", "50216", "50516"]
# folders = ["46916", "47216"]
nSpace = 12
inside = nSpace // 2 - 1

outv = 12 * 11 * 11
inv = 828

outmaxomegalst = []
maxomegalst = []
maxvlst = []
suspOmealllst = []
suspValllst = []
xrecorded = False

for f in folders:
    polyalstin = []
    polyalstine = []
    susplstin = []
    susplstine = []
    polyalstout = []
    polyalstoute = []
    susplstout = []
    susplstoute = []
    # start
    fr = filenames.format(f, f, "R")
    pr = np.real(LoadMathematicaCSV(fr))
    nEdge = nSpace // 2 - 1
    prSq = np.round(pr * pr)
    numSite = np.array([0 for _ in range(len(pr))])
    for x in range(-nEdge, nEdge + 1):
        for y in range(-nEdge, nEdge + 1):
            idxInPrSq = prSq.tolist().index(x * x + y * y)
            numSite[idxInPrSq] = numSite[idxInPrSq] + 1
    fin = filenames.format(f, f, "In")
    fout = filenames.format(f, f, "Out")
    # r-data
    xlst = []
    plst = []
    pelst = []
    slst = []
    selst = []
    for ome in range(21):
        # in
        din = LoadMathematicaCSV(fin)
        arrin = np.real(din[ome, :]) if useRealInsteadOfAbs else np.abs(din[ome, :])
        # if 0 == ome:
        #     print(f + ": ", np.shape(arrin))
        v, s = JacknifeMean(arrin)
        _, _, t = AutoCorrelationSingleVariable(arrin)
        polyalstin.append(v)
        polyalstine.append(s * np.sqrt(2 * t))
        v, s = JacknifeCumulant(arrin)
        susplstin.append(v)
        susplstine.append(s * np.sqrt(2 * t))
        # out
        dout = LoadMathematicaCSV(fout)
        arrout = np.real(dout[ome, :]) if useRealInsteadOfAbs else np.abs(dout[ome, :])
        v, s = JacknifeMean(arrout)
        _, _, t = AutoCorrelationSingleVariable(arrout)
        polyalstout.append(v)
        polyalstoute.append(s * np.sqrt(2 * t))
        v, s = JacknifeCumulant(arrout)
        susplstout.append(v)
        susplstoute.append(s * np.sqrt(2 * t))
        # r
        fo = filenames.format(f, f, "O{}".format(ome))
        dome = np.real(LoadMathematicaCSV(fo)) if useRealInsteadOfAbs else np.abs(LoadMathematicaCSV(fo))
        # print("working for ome {}".format(ome))
        for r in range(len(pr)):
            if 0 == r:
                continue
            if pr[r] >= inside:
                continue
            ro = ome * pr[r]
            polyad = dome[:, r]
            v, s = JacknifeMean(polyad)
            _, _, t = AutoCorrelationSingleVariable(polyad)
            xlst.append(ro)
            plst.append(v)
            pelst.append(s * np.sqrt(2 * t))
            v, s = JacknifeCumulant(polyad)
            slst.append(numSite[r] * v)
            selst.append(numSite[r] * s * np.sqrt(2 * t))
    polyalstin = np.array(polyalstin) / 3
    # print(polyalstin)
    polyalstine = np.array(polyalstine) / 3
    # print(polyalstine)
    susplstin = inv * np.array(susplstin) / 9
    # print(susplstin)
    susplstine = inv * np.array(susplstine) / 9
    # print(susplstine)
    polyalstout = np.array(polyalstout) / 3
    # print(polyalstout)
    polyalstoute = np.array(polyalstoute) / 3
    # print(polyalstoute)
    susplstout = outv * np.array(susplstout) / 9
    # print(susplstout)
    susplstoute = outv * np.array(susplstoute) / 9
    # print(susplstoute)
    xlst = np.array(xlst)
    plst = np.array(plst)
    pelst = np.array(pelst)
    slst = np.array(slst) * 4 / 3 # it is slst * Nz / 9 with Nz=12, Nt=3
    selst = np.array(selst) * 4 / 3
    outmaxomegalst.append(np.argmax(susplstout))
    maxomegalst.append(np.argmax(susplstin))
    maxvlst.append(np.argmax(slst))
    suspOmealllst.append(susplstin)
    suspOmealllst.append(susplstine)
    suspOmealllst.append(susplstout)
    suspOmealllst.append(susplstoute)
    if not xrecorded:
        xrecorded = True
        suspValllst.append(xlst)
    suspValllst.append(slst)
    suspValllst.append(selst)
    print(xlst[np.argmax(slst)])
    # print(xlst)
    # print(plst)
    # print(pelst)
    # print(slst)
    # print(selst)
    errorbar(range(len(susplstin)), [susplstin / 9, susplstout / 9], [susplstine / 9, susplstoute / 9])
    # errorbar(xlst, [plst / 3], [pelst / 3], markers=[MarkerStyle.circle], linestyles=[LineStyle.none])
    errorbar(xlst, [slst], [selst], markers=[MarkerStyle.circle], linestyles=[LineStyle.none])

outmaxomegalst = np.reshape(np.array(outmaxomegalst), (len(folders), 1))
maxomegalst = np.reshape(np.array(maxomegalst), (len(folders), 1))
maxvlst = np.reshape(np.array(maxvlst), (len(folders), 1))
np.savetxt(writefolder + "maxomega.csv", np.hstack((outmaxomegalst, maxomegalst, maxvlst)), fmt='%d', delimiter=',')
np.savetxt(writefolder + "susp.csv", suspOmealllst, delimiter=',')
np.savetxt(writefolder + "suspv.csv", suspValllst, delimiter=',')
