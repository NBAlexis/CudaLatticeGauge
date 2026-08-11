from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

hd = "580"
folderHeads = f"H:\\MesonMassElectric\\{hd}01\\Chiral\\"
# folderHeads = "H:\\BS01Polyakov\\"

cvu = []
svu = []
tvu = []
cvd = []
svd = []
tvd = []

for i in range(8):
    fileNames = folderHeads + f"EC{hd}_01_M0__{i}_condensatepCCHeavyChiralKS.csv"
    testarray = LoadMathematicaCSV(fileNames)
    testarrayre = np.real(testarray[100:])
    print(np.shape(testarrayre))
    _, _, t = AutoCorrelationSingleVariable(testarrayre, s=2)
    v, s = JacknifeMean(testarrayre)
    cvu.append(-v)
    svu.append(s * np.sqrt(2 * t))
    tvu.append(t)
    fileNames = folderHeads + f"EC{hd}_01_M0__{i}_condensatepCCLightChiralKS.csv"
    testarray = LoadMathematicaCSV(fileNames)
    testarrayre = np.real(testarray[100:])
    _, _, t = AutoCorrelationSingleVariable(testarrayre, s=2)
    v, s = JacknifeMean(testarrayre)
    cvd.append(-v)
    svd.append(s * np.sqrt(2 * t))
    tvd.append(t)


print(PrintAsMatlabArray(cvu, f"cvd{hd}"))
print(PrintAsMatlabArray(svu, f"svd{hd}"))
print(PrintAsMatlabArray(tvu, f"tvd{hd}"))

print(PrintAsMatlabArray(cvd, f"cvu{hd}"))
print(PrintAsMatlabArray(svd, f"svu{hd}"))
print(PrintAsMatlabArray(tvd, f"tvu{hd}"))

"""
cvd530=[0.42627253,0.42660263,0.43009567,0.43291419,0.436655,0.44088493,0.44608765,0.45180452]
svd530=[0.00052638,0.00029883,0.00037645,0.00044034,0.00044269,0.00031756,0.00035169,0.00038535]
tvd530=[19.79237282,7.84197827,12.3086356,15.09677172,14.66353677,10.10362549,11.38249014,15.88974741]
cvu530=[0.42627132,0.42761764,0.43370347,0.44068277,0.44932995,0.45929342,0.47017013,0.48155151]
svu530=[0.00050581,0.0003333,0.0003548,0.00039923,0.00043169,0.00029852,0.00030248,0.00034255]
tvu530=[18.61248365,9.29181979,11.31942966,12.65581393,14.98259657,9.86354541,10.1627979,14.28715974]

cvd580=[0.12347507,0.13738328,0.14096185,0.14580295,0.15145409,0.15777233,0.1643219,0.17110026]
svd580=[6.98143352e-05,8.31617459e-05,8.80003112e-05,7.41925414e-05,7.21382440e-05,9.94330171e-05,1.43751073e-04,1.06963350e-04]
tvd580=[6.79237435,7.70775779,8.4487072,5.57329702,4.62774889,8.09441909,14.90785518,7.96171524]
cvu580=[0.12347391,0.14092566,0.1512832,0.16367342,0.17743097,0.19215992,0.20713704,0.22174896]
svu580=[6.74151709e-05,7.43958069e-05,8.19806056e-05,7.21871604e-05,8.10966879e-05,1.05368545e-04,1.51704383e-04,1.22887485e-04]
tvu580=[6.26399024,5.87395242,6.42512971,4.55097561,4.41869839,6.47278167,11.19868449,6.65100603]
"""
