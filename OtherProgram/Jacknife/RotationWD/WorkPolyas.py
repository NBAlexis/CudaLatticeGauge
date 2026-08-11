import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMatlabArray
from Visualization import errorbar

"""
counts:
654
1544
3050
5244
"""

counts = [654, 13*13*6, 1544, 17*17*8, 3050, 21*21*10, 5244, 25*25*12]

f1 = ["H:\\RotationWD\\25153Nt3\\Polyakov\\NC25153__polyakov_Nt3_In.csv"]
f2 = ["H:\\RotationWD\\25153Nt3\\Polyakov\\NC25153__polyakov_Nt3_Out.csv"]
f3 = ["H:\\RotationWD\\502161Nt4\\Polyakov\\NC502161__polyakov_Nt4_In.csv"]
f4 = ["H:\\RotationWD\\502161Nt4\\Polyakov\\NC502161__polyakov_Nt4_Out.csv"]
f5 = ["H:\\RotationWD\\5283161Nt5\\Polyakov\\NC5283161__polyakov_Nt5_In.csv"]
f6 = ["H:\\RotationWD\\5283161Nt5\\Polyakov\\NC5283161__polyakov_Nt5_Out.csv"]
f7 = ["H:\\RotationWD\\535161Nt6\\NewPolyakov\\NC535161__polyakov_Nt6_In.csv"]
f8 = ["H:\\RotationWD\\535161Nt6\\NewPolyakov\\NC535161__polyakov_Nt6_Out.csv"]

datain = None
useRealOrAbs = True
for file in f8:
    if datain is None:
        datain = LoadMathematicaCSV(file)
    else:
        datain = np.hstack((datain, LoadMathematicaCSV(file)))
print(np.shape(datain))

polyalst = []
polyalste = []
susplst = []
susplste = []

for i in range(0, 11):
    arrin = np.real(datain[i, :]) if useRealOrAbs else np.abs(datain[i, :])
    v, s = JacknifeMean(arrin)
    _, _, t = AutoCorrelationSingleVariable(arrin)
    polyalst.append(v)
    polyalste.append(s * np.sqrt(2 * t))
    v, s = JacknifeCumulant(arrin)
    susplst.append(v)
    susplste.append(s * np.sqrt(2 * t))

polyalst = np.array(polyalst)
polyalste = np.array(polyalste)
susplst = np.array(susplst)
susplste = np.array(susplste)

print(PrintAsMatlabArray(polyalst / 3, "p"))
print(PrintAsMatlabArray(polyalste / 3, "pe"))
print(PrintAsMatlabArray(susplst / 9, "susp"))
print(PrintAsMatlabArray(susplste / 9, "suspe"))


errorbar(range(len(polyalst)), [polyalst / 3], [polyalste / 3])
errorbar(range(len(susplst)), [susplst / 9], [susplste / 9])

"""
count14=0
for xx in range(13):
    x = xx - 6
    for yy in range(13):
        y = yy - 6
        if x*x + y*y < 36:
            count14=count14+1
count14=count14*6
print(count14)

count14=0
for xx in range(17):
    x = xx - 8
    for yy in range(17):
        y = yy - 8
        if x*x + y*y < 64:
            count14=count14+1
count14=count14*8
print(count14)

count14=0
for xx in range(21):
    x = xx - 10
    for yy in range(21):
        y = yy - 10
        if x*x + y*y < 100:
            count14=count14+1
count14=count14*10
print(count14)

count14=0
for xx in range(25):
    x = xx - 12
    for yy in range(25):
        y = yy - 12
        if x*x + y*y < 144:
            count14=count14+1
count14=count14*12
print(count14)
"""
