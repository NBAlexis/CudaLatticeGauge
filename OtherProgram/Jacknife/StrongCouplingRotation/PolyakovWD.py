import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray, \
    PrintAsMatlabArray
from Visualization import errorbar

paramlst = ["0521", "10209", "2020"]

header = "I:\\StrongCoupling\\WilsonDirac\\Polyakov\\"

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

for paramidx in range(len(paramlst)):
    thisxiv = []
    thisxis = []
    thisxit = []
    thisxive = []
    thisxise = []
    lastone = []
    for omega in range(11):
        data = np.abs(LoadMathematicaCSV(f"{header}NC{paramlst[paramidx]}__{omega}_polyakov_In.csv"))
        lastone.append(data[len(data) - 1])
        v, s = JacknifeMean(data)
        _, _, t = AutoCorrelationSingleVariable(data)
        thisxiv.append(v)
        thisxit.append(t)
        thisxive.append(s * np.sqrt(2 * t))
        v, s = JacknifeCumulant(data, f"susp{paramlst[paramidx]}-{omega}")
        thisxis.append(v)
        thisxise.append(s * np.sqrt(2 * t))
    errorbar(range(len(thisxiv)), [thisxiv], [thisxive])
    allv.append(thisxiv)
    allve.append(thisxive)
    alls.append(thisxis)
    allse.append(thisxise)
    allt.append(thisxit)
    thisxiv = []
    thisxis = []
    thisxit = []
    thisxive = []
    thisxise = []
    lastone = []
    for omega in range(11):
        data = np.abs(LoadMathematicaCSV(f"{header}NC{paramlst[paramidx]}__{omega}_polyakov.csv"))
        lastone.append(data[len(data) - 1])
        v, s = JacknifeMean(data)
        _, _, t = AutoCorrelationSingleVariable(data)
        thisxiv.append(v)
        thisxit.append(t)
        thisxive.append(s * np.sqrt(2 * t))
        v, s = JacknifeCumulant(data, f"susp{paramlst[paramidx]}-{omega}")
        thisxis.append(v)
        thisxise.append(s * np.sqrt(2 * t))
    print(PrintAsMatlabArray(lastone, "lastone"))
    allvo.append(thisxiv)
    allveo.append(thisxive)
    allso.append(thisxis)
    allseo.append(thisxise)
    allto.append(thisxit)

allv = np.array(allv) / 3
allve = np.array(allve) / 3
alls = np.array(alls) / 9
allse = np.array(allse) / 9
allvo = np.array(allvo) / 3
allveo = np.array(allveo) / 3
allso = np.array(allso) / 9
allseo = np.array(allseo) / 9

alltosave = np.vstack((allv, allve, alls, allse, np.array(allt), allvo, allveo, allso, allseo, np.array(allto)))
np.savetxt(f"Data/PolyaWD.csv", alltosave, delimiter=",")

