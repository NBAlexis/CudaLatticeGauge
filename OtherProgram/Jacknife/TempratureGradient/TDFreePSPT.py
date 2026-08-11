import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray
from Visualization import errorbar

header = "I:\\BetaGradient\\TemperatureDist\\TDFree\\Wilson\\"
nt = 6
# betalst = ["565", "570", "575", "580", "585", "590", "595", "600", "605", "610", "615"]

# we mistake store all results together
testarray = np.real(LoadMathematicaCSV(header + "TDFree6__615_wilsonloops.csv"))

valuePS = []
stdPS = []
valuesPS = []
stdsPS = []
tsPS = []

valuePT = []
stdPT = []
valuesPT = []
stdsPT = []
tsPT = []

for bt in range(11):
    ress = testarray[bt * 2900: (bt + 1) * 2900, [0, 1, 3]]
    rest = testarray[bt * 2900: (bt + 1) * 2900, [2, 4, 5]]
    ress = np.sum(ress, axis=1) / 3
    rest = np.sum(rest, axis=1) / 3
    print(np.shape(ress))
    v, s = JacknifeMean(ress)
    _, _, t = AutoCorrelationSingleVariable(ress)
    valuePS.append(v)
    tsPS.append(t)
    stdPS.append(s* np.sqrt(2 * t))
    v, s = JacknifeCumulant(ress)
    valuesPS.append(v)
    stdsPS.append(s* np.sqrt(2 * t))
    # PT
    v, s = JacknifeMean(rest)
    _, _, t = AutoCorrelationSingleVariable(rest)
    valuePT.append(v)
    tsPT.append(t)
    stdPT.append(s* np.sqrt(2 * t))
    v, s = JacknifeCumulant(rest)
    valuesPT.append(v)
    stdsPT.append(s* np.sqrt(2 * t))

valuePS = np.array(valuePS)
print(PrintAsMathematicaArray(3 - valuePS, "PS"))
print(PrintAsMathematicaArray(stdPS, "PSe"))
print(PrintAsMathematicaArray(valuesPS, "chiPS"))
print(PrintAsMathematicaArray(stdsPS, "chiPSe"))
print(PrintAsMathematicaArray(tsPS, "PSt"))

valuePT = np.array(valuePT)
print(PrintAsMathematicaArray(3 - valuePT, "PT"))
print(PrintAsMathematicaArray(stdPT, "PTe"))
print(PrintAsMathematicaArray(valuesPT, "chiPT"))
print(PrintAsMathematicaArray(stdsPT, "chiPTe"))
print(PrintAsMathematicaArray(tsPT, "PTt"))

"""

PS={1.22693186,1.20259272,1.18131493,1.16212281,1.14462605,1.1284001,1.11329309,1.09887944,1.08522268,1.07227985,1.05968839};
PSe={1.48742508*^-04,1.14210347*^-04,1.20181158*^-04,1.13467584*^-04,8.95261608*^-05,1.05173773*^-04,7.68594537*^-05,8.35553787*^-05,7.31608979*^-05,7.90547666*^-05,8.13592151*^-05};
chiPS={4.82647303*^-06,4.21603859*^-06,4.04028930*^-06,3.83093354*^-06,3.53755872*^-06,3.47034826*^-06,3.04541822*^-06,3.29704584*^-06,2.93635365*^-06,2.75628921*^-06,2.94574207*^-06};
chiPSe={4.47291081*^-07,3.17232033*^-07,3.42093419*^-07,3.14978397*^-07,2.33771891*^-07,2.79508582*^-07,1.79193075*^-07,2.11482017*^-07,1.82209693*^-07,1.82773018*^-07,2.01612237*^-07};
PSt={6.6467321,4.48615569,5.18356421,4.87311605,3.28521985,4.62178696,2.81264941,3.07037724,2.64312271,3.28775726,3.25826646};
PT={1.42655359,1.39500184,1.36838677,1.34458157,1.32298205,1.3028312,1.28441886,1.26727332,1.25115974,1.23562407,1.22103196};
PTe={1.92397946*^-04,1.65633382*^-04,1.54520902*^-04,1.07320497*^-04,1.32738759*^-04,1.30525112*^-04,9.99249887*^-05,9.09735948*^-05,9.38285053*^-05,9.33179361*^-05,9.94663252*^-05};
chiPT={7.52438135*^-06,6.60919449*^-06,5.96330556*^-06,5.50813434*^-06,5.22914517*^-06,5.19316086*^-06,4.53479749*^-06,4.30421554*^-06,4.21829863*^-06,4.16314367*^-06,3.95204240*^-06};
chiPTe={7.25019854*^-07,5.89675220*^-07,4.99508633*^-07,3.53894105*^-07,4.10237544*^-07,3.90601897*^-07,3.00423019*^-07,2.72347240*^-07,2.78549715*^-07,2.76729000*^-07,2.86011383*^-07};
PTt={7.1334249,6.01887404,5.80571328,3.03199849,4.88576943,4.75690364,3.19270277,2.78807807,3.02621651,3.03303347,3.62993147};

This result is bad to explain the `xi' values...

"""