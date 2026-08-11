import numpy as np

from JacknifePrograms import LoadMathematicaCSV, FitWDMass, JacknifeWilsonLoop, PrintAsMathematicaArray

f1 = "H:\\50516Matching\\NC50516_pion.csv"
f2 = "H:\\50516Matching\\NC50516_rho.csv"
# f3 = "F:\\Builds\\DoubleBuild-07-25\\Release\\NC535161__VR_R.csv"
# f4 = "F:\\Builds\\DoubleBuild-07-25\\Release\\NC535161__VR.csv"
data1 = np.real(LoadMathematicaCSV(f1))
print(FitWDMass(data1))
data2 = np.real(LoadMathematicaCSV(f2))
print(FitWDMass(data2))
# data3 = np.real(LoadMathematicaCSV(f3))
# data4 = np.real(LoadMathematicaCSV(f4))
# r0, r0e, r1, r1e, c0, c0e = JacknifeWilsonLoop(data4, data3, 24, 1, 3, 10, "", True)
# print(PrintAsMathematicaArray([r0, r0e, r1, r1e, c0, c0e]))