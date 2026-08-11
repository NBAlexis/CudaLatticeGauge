import numpy as np

from JacknifePrograms import PrintAsMathematicaArray

f530 = np.loadtxt("../DrawData/53001_effmass.csv", delimiter=",")

print(PrintAsMathematicaArray(f530))