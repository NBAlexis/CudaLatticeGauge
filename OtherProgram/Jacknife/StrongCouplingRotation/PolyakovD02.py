import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray, \
    PrintAsMatlabArray
from Visualization import errorbar

# xilst = [14, 16, 18, 20, 21, 22, 23, 24, 242, 245, 25, 26, 28, 30, 32, 34]
# xilst = [226, 228, 230, 232, 237, 240, 243, 246]
xilst = [2270, 2275, 2280, 2285, 2290, 2295, 2300, 2305, 2310, 2315, 2320, 2325, 2330, 2335, 2340]
loadfolder = "QP02"
fileheader = "QP020"
csvtosave = "QP02.csv"
# xilst = [18, 22, 34, 38]
# xilst = [25, 30, 35, 37, 40, 45, 50]
storesusp = True
header = f"I:\\StrongCoupling\\{loadfolder}\\Polyakov\\"
# header = "F:\\Cppwork\\CudaLatticeGauge\\Bin\\Release\\QP02\\Polyakov\\"
maxomega = 31

allv = []
allve = []
alls = []
allse = []
allt = []

allvo = []
allveo = []
allso = []
allseo = []
allto = []

for xiidx in range(len(xilst)):
    thisxiv = []
    thisxis = []
    thisxit = []
    thisxive = []
    thisxise = []
    for omega in range(maxomega):
        data = np.abs(LoadMathematicaCSV(f"{header}{fileheader}{xilst[xiidx]}__{omega}_polyakov_In.csv"))
        if storesusp:
            v, s = JacknifeMean(data)
            _, _, t = AutoCorrelationSingleVariable(data)
            thisxiv.append(v)
            thisxit.append(t)
            thisxive.append(s * np.sqrt(2 * t))
            v, s = JacknifeCumulant(data, f"susp{xilst[xiidx]}-{omega}")
            thisxis.append(v)
            thisxise.append(s * np.sqrt(2 * t))
        else:
            v, ve, t = AutoCorrelationSingleVariable(data)
            print(f"polya{xilst[xiidx]}-{omega}")
            thisxiv.append(v)
            thisxit.append(t)
            thisxive.append(ve)
    errorbar(range(len(thisxiv)), [thisxiv], [thisxive])
    if storesusp:
        errorbar(range(len(thisxis)), [thisxis], [thisxise])
    allv.append(thisxiv)
    allve.append(thisxive)
    if storesusp:
        alls.append(thisxis)
        allse.append(thisxise)
    allt.append(thisxit)
    thisxiv = []
    thisxis = []
    thisxit = []
    thisxive = []
    thisxise = []
    lastone = []
    for omega in range(maxomega):
        data = np.abs(LoadMathematicaCSV(f"{header}{fileheader}{xilst[xiidx]}__{omega}_polyakov.csv"))
        lastone.append(data[len(data) - 1])
        if storesusp:
            v, s = JacknifeMean(data)
            _, _, t = AutoCorrelationSingleVariable(data)
            thisxiv.append(v)
            thisxit.append(t)
            thisxive.append(s * np.sqrt(2 * t))
            v, s = JacknifeCumulant(data, f"susp{xilst[xiidx]}-{omega}")
            thisxis.append(v)
            thisxise.append(s * np.sqrt(2 * t))
        else:
            v, ve, t = AutoCorrelationSingleVariable(data)
            print(f"polya{xilst[xiidx]}-{omega}")
            thisxiv.append(v)
            thisxit.append(t)
            thisxive.append(ve)
    print(PrintAsMatlabArray(lastone, "lastone"))
    allvo.append(thisxiv)
    allveo.append(thisxive)
    if storesusp:
        allso.append(thisxis)
        allseo.append(thisxise)
    allto.append(thisxit)

allv = np.array(allv) / 3
allve = np.array(allve) / 3
if storesusp:
    alls = np.array(alls) / 9
    allse = np.array(allse) / 9
allvo = np.array(allvo) / 3
allveo = np.array(allveo) / 3
if storesusp:
    allso = np.array(allso) / 9
    allseo = np.array(allseo) / 9

if storesusp:
    alltosave = np.vstack((allv, allve, alls, allse, np.array(allt), allvo, allveo, allso, allseo, np.array(allto)))
else:
    alltosave = np.vstack((allv, allve, np.array(allt), allvo, allveo, np.array(allto)))
np.savetxt(f"Data/{csvtosave}", alltosave, delimiter=",")

"""

"""