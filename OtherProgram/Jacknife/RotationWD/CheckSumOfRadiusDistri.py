import numpy as np

from JacknifePrograms import LoadMathematicaCSV
from RotationWD.UsefulFunctions import getPrAndNumberOfSites

diskName = "H:\\"

omega = 8
checkf1 = diskName + "RotationWD\\5283161Nt5\\NewChiral\\NC5283161__condensateChiral_Nt5_O{}.csv".format(omega)
checkf2 = diskName + "RotationWD\\5283161Nt5\\NewChiral\\NC5283161__condensateChiral_Nt5_All_O{}.csv".format(omega)

_, numpr, allsitenum, _ = getPrAndNumberOfSites(22)


d1 = LoadMathematicaCSV(checkf1)
d1 = np.dot(d1, numpr) / allsitenum
print(d1)
d2 = LoadMathematicaCSV(checkf2)
print(d2)
print(np.sum(np.abs(d1 - d2)))