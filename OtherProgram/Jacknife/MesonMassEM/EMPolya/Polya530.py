from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

hd = "580"
folderHeads = f"H:\\MesonMassElectric\\{hd}01\\Polyakov\\"

pv = []
ps = []
pt = []
for i in range(8):
    fileNames = folderHeads + f"EC{hd}_01_M0__{i}_polyakov.csv"
    testarray = LoadMathematicaCSV(fileNames)
    testarray = np.abs(testarray[100:]) / 3
    _, _, t = AutoCorrelationSingleVariable(testarray)
    v, s = JacknifeMean(testarray)
    pv.append(v)
    ps.append(s * np.sqrt(2 * t))
    pt.append(t)

print(PrintAsMatlabArray(pv, f"pv{hd}"))
print(PrintAsMatlabArray(ps, f"ps{hd}"))
print(PrintAsMatlabArray(pt, f"pt{hd}"))

"""
pv530=[0.0033201,0.0026912,0.00270117,0.00276212,0.00267603,0.00275637,0.00274342,0.00274035]
ps530=[5.51044404e-05,3.90601580e-05,4.10175102e-05,4.04781424e-05,4.13719362e-05,3.89722533e-05,3.99380658e-05,3.89250503e-05]
pt530=[0.68916761,0.51998362,0.57504206,0.51511144,0.57557797,0.51866411,0.52623795,0.5]


pv580=[0.11883079,0.09567673,0.09877896,0.10036176,0.0983406,0.09654088,0.09502173,0.09652667]
ps580=[0.00051487,0.00144262,0.00091483,0.00061934,0.00047942,0.00100016,0.00081058,0.0007859]
pt580=[10.81493853,47.81381496,20.70256406,13.06086484,8.85416287,28.49871172,19.10098587,16.52552985]

"""


"""
cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(8):
    fileNames = folderHeads + f"EC530_01_M0__{i}_polyakov_XSlice.csv"
    testarray = LoadMathematicaCSV(fileNames)
    testarrayre = np.real(testarray)
    testarrayim = np.imag(testarray)
    resre = []
    resim = []
    for j in range(24):
        v, _ = JacknifeMean(testarrayre[:,j])
        resre.append(v)
        v, _ = JacknifeMean(testarrayim[:,j])
        resim.append(v)
    print(resre)
    print(resim)
    plt.plot(resre, resim)
    plt.show()
"""

# print(PrintAsMathematicaArray(values, "susp"))
# print(PrintAsMathematicaArray(stds, "suspe"))
# print(PrintAsMathematicaArray(tvalues, "t"))

# errorbar(range(len(values)), [cv], [sv])
# errorbar(range(len(values)), [values], [stds])

