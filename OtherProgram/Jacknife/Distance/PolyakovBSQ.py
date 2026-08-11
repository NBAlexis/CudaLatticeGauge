from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads = "H:\\Quench242424-6\\Polyakov\\"


lst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]

cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(len(lst)):
    fileNames = folderHeads + "BSQ__{}_polyakov.csv".format(lst[i])
    testarray = LoadMathematicaCSV(fileNames)
    testarrayre = np.abs(testarray)
    v, s = JacknifeMean(testarrayre, lst[i])
    _, _, t = AutoCorrelationSingleVariable(testarrayre)
    cv.append(v / 3)
    sv.append(s * np.sqrt(2 * t) / 3)
    v, s = JacknifeCumulant(testarrayre, lst[i])
    print("==============", i, "with {} results".format(len(testarrayre)))
    print(v)
    print(s)
    values.append(v / 9)
    stds.append(s * np.sqrt(2 * t) / 9)
    tvalues.append(t)

print(PrintAsMathematicaArray(cv, "v"))
print(PrintAsMathematicaArray(sv, "s"))
print(PrintAsMathematicaArray(tvalues, "t"))
print(PrintAsMathematicaArray(values, "susp"))
print(PrintAsMathematicaArray(stds, "suspe"))


errorbar(range(len(values)), [cv], [sv])
errorbar(range(len(values)), [values], [stds])

"""
v={0.00340909,0.00389511,0.00480925,0.00595234,0.01011893,0.0366135,0.08127074,0.09399013,0.10454116,0.11333382,0.12333582};
s={6.00490349*^-05,8.67482800*^-05,2.03772431*^-04,3.37164912*^-04,1.27473983*^-03,8.25899939*^-03,2.68612509*^-03,9.12186119*^-04,6.76750843*^-04,5.34846718*^-04,4.70770318*^-04};
t={1.58765448,2.69872513,9.17451221,17.53497568,85.05886816,343.4728837,132.49585696,45.14079142,26.33486194,17.64166468,14.40242813};
susp={3.29324524*^-06,4.04325465*^-06,6.56259909*^-06,9.40042698*^-06,2.77007493*^-05,2.87958839*^-04,7.89620058*^-05,2.67279563*^-05,2.52170667*^-05,2.35118666*^-05,2.23126129*^-05};
suspe={1.62982598*^-07,2.49277020*^-07,7.94441840*^-07,1.51183953*^-06,9.73144412*^-06,1.27159686*^-04,3.82292351*^-05,6.71062231*^-06,4.73854379*^-06,3.53514467*^-06,3.15259781*^-06};
"""