import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray
from Visualization import errorbar

header = "I:\\BetaGradient\\TemperatureDist\\TDFree\\Polyakov\\"
nt = 6
betalst = ["565", "570", "575", "580", "585", "590", "595", "600", "605", "610", "615"]

value = []
std = []
values = []
stds = []
ts = []
for bt in betalst:
    filename = header + f'TDFree6__{bt}_polyakov.csv'
    testarray = np.abs(LoadMathematicaCSV(filename))
    v, s = JacknifeMean(testarray, bt)
    _, _, t = AutoCorrelationSingleVariable(testarray)
    value.append(v / 3)
    ts.append(t)
    std.append(s* np.sqrt(2 * t) / 3)
    v, s = JacknifeCumulant(testarray)
    values.append(v / 9)
    stds.append(s* np.sqrt(2 * t) / 9)


print(PrintAsMathematicaArray(value))
print(PrintAsMathematicaArray(std))
print(PrintAsMathematicaArray(values))
print(PrintAsMathematicaArray(stds))
print(PrintAsMathematicaArray(ts))
errorbar(range(len(values)), [value], [std])
errorbar(range(len(values)), [values], [stds])

"""

{0.00354348,0.00417738,0.00598724,0.00984061,0.02580634,0.04856973,0.0664534,0.07560518,0.08391144,0.09195023,0.09890047};
{5.91756120*^-05,9.86642410*^-05,5.59883057*^-04,2.29194942*^-03,5.61295916*^-03,5.15551663*^-03,1.32900723*^-03,6.91126790*^-04,1.01513220*^-03,8.13845321*^-04,8.06614268*^-04};
{3.35769562*^-06,4.40803983*^-06,9.52537722*^-06,3.00625745*^-05,1.48284042*^-04,1.17758332*^-04,3.16934558*^-05,2.28339720*^-05,2.66885836*^-05,2.57940555*^-05,2.87771830*^-05};
{1.56344414*^-07,2.97543872*^-07,2.29433991*^-06,1.39986739*^-05,7.90817329*^-05,6.47872566*^-05,9.85461399*^-06,4.75985016*^-06,7.28344997*^-06,6.66798949*^-06,5.86262153*^-06};
{1.51221031,3.2021528,47.71780629,253.36807451,308.07563523,327.28096142,80.80776413,30.33206616,55.9870629,37.23335018,32.78321361};

"""