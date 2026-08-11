import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMatlabArray
from Visualization import errorbar, LineStyle, MarkerStyle

diskName = "H:\\"
f1 = diskName + "RotationWD\\25153Nt3\\Polyakov\\NC25153__polyakov_Nt3_In.csv"
f2 = diskName + "RotationWD\\25153Nt3\\Polyakov\\NC25153__polyakov_Nt3_Out.csv"
f3 = diskName + "RotationWD\\5283161Nt5\\Polyakov\\NC5283161__polyakov_Nt5_{}.csv"
f4 = diskName + "RotationWD\\502161Nt4\\Polyakov\\NC502161__polyakov_Nt4_Out.csv"

filename = f3
rdata = np.real(LoadMathematicaCSV(filename.format("R")))
pdata = np.real(LoadMathematicaCSV(filename.format("O0")))
print(np.shape(pdata))
prlst = []
prelst = []
for i in range(len(rdata)):
    pr = pdata[:, i]
    v, s = JacknifeMean(pr)
    _, _, t = AutoCorrelationSingleVariable(pr)
    prlst.append(v)
    prelst.append(s * np.sqrt(2 * t))

prlst = np.array(prlst) / 3
prelst = np.array(prelst) / 3

errorbar(rdata, [prlst], [prelst], linestyles=[LineStyle.none], markers=[MarkerStyle.circle])

print(PrintAsMatlabArray(rdata))
print(PrintAsMatlabArray(prlst))
print(PrintAsMatlabArray(prelst))

