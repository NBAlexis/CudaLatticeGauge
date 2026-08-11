import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray
from Visualization import errorbar

header = "I:\\BetaGradient\\TemperatureDist\\TDFree\\Boson\\"
nt = 6
betalst = ["565", "570", "575", "580", "585", "590", "595", "600", "605", "610", "615"]

value = []
std = []
values = []
stds = []
ts = []
for bt in betalst:
    filename = header + f'TDFree{nt}__{bt}_bosonvalueR.csv'
    testarray = np.real(LoadMathematicaCSV(filename))
    v, s = JacknifeMean(testarray, bt)
    _, _, t = AutoCorrelationSingleVariable(testarray)
    value.append(v)
    ts.append(t)
    std.append(s* np.sqrt(2 * t))
    v, s = JacknifeCumulant(testarray)
    values.append(v)
    stds.append(s* np.sqrt(2 * t))

print(PrintAsMathematicaArray(value))
print(PrintAsMathematicaArray(std))
print(PrintAsMathematicaArray(values))
print(PrintAsMathematicaArray(stds))
print(PrintAsMathematicaArray(ts))
errorbar(range(len(values)), [value], [std])
errorbar(range(len(values)), [values], [stds])

"""

{1.14538177,1.14468307,1.14441882,1.14416918,1.14398712,1.14369718,1.14358023,1.14359298,1.14366219,1.14362215,1.14378343};
{6.07238964*^-05,6.22013825*^-05,6.51477294*^-05,7.48639875*^-05,6.73160991*^-05,6.87217318*^-05,5.94263161*^-05,6.58680413*^-05,6.27419042*^-05,6.42164694*^-05,7.13843243*^-05};
{2.39561570*^-06,2.38850717*^-06,2.44299726*^-06,2.54328609*^-06,2.54506750*^-06,2.52870239*^-06,2.27331788*^-06,2.41455695*^-06,2.39773448*^-06,2.38067236*^-06,2.48055510*^-06};
{1.30411445*^-07,1.35823653*^-07,1.38720832*^-07,1.78094486*^-07,1.45522557*^-07,1.52264170*^-07,1.23858246*^-07,1.50360496*^-07,1.38871485*^-07,1.45226401*^-07,1.66827285*^-07};
{2.23187637,2.34877489,2.51908895,3.19535263,2.58170332,2.70806,2.25250288,2.60543311,2.38057784,2.5116619,2.97868657};


"""