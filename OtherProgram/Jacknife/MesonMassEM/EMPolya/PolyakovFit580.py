from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

from MesonMassEM.Mass.CorrelatedFitGPT import PolyaFit, polyakov_fit, PolyaFit2, polyakov_fit2, polyakov_fit580b, \
    PolyaFit580, PolyaFit5802, polyakov_fit580b2

hd = "580"
folderHeads = f"H:\\MesonMassElectric\\{hd}01\\Polyakov\\"

pt580 = [5.024028375384501, 41.126456042762996, 27.06598166119774, 10.08989122533077, 7.944354294370299, 10.121061163852499, 16.983199382378825, 10.954232229869973]

cv = []
sv = []
values = []
stds = []
tvalues = []


fileNames = folderHeads + f"EC{hd}_01_M0__1_polyakov_XSlice.csv"
testarray = LoadMathematicaCSV(fileNames)[100:, :]
a, ae, b, be, c, ce, d, de, e, ee, chi2 = PolyaFit580(testarray)
print(a, ae * np.sqrt(2 * pt580[1]), b, be * np.sqrt(2 * pt580[1]), c, ce * np.sqrt(2 * pt580[1]), d, de * np.sqrt(2 * pt580[1]), e, ee * np.sqrt(2 * pt580[1]), chi2 / (2 * pt580[1]))
testarrayre = np.real(testarray)
testarrayim = np.imag(testarray)
resre = []
resres = []
resim = []
resims = []
for j in range(24):
    _, _, t = AutoCorrelationSingleVariable(testarrayre[:, j], s=2.0)
    v, s = JacknifeMean(testarrayre[:, j])
    resre.append(v)
    resres.append(s * np.sqrt(2 * t))
    _, _, t = AutoCorrelationSingleVariable(testarrayim[:, j], s=2.0)
    v, s = JacknifeMean(testarrayim[:, j])
    resim.append(v)
    resims.append(s * np.sqrt(2 * t))
print(PrintAsMatlabArray(resre))
print(PrintAsMatlabArray(resres))
print(PrintAsMatlabArray(resim))
print(PrintAsMatlabArray(resims))
plt.plot(resre, resim)
xx = np.array([m * 0.1 for m in range(240)])
yy = polyakov_fit580b(xx, [a, b, c, d, e])
yyre = np.real(yy)
yyim = np.imag(yy)
plt.plot(yyre, yyim)
plt.show()

fileNames = folderHeads + f"EC{hd}_01_M0__2_polyakov_XSlice.csv"
testarray = LoadMathematicaCSV(fileNames)
a, ae, b, be, c, ce, d, de, e, ee, chi2 = PolyaFit5802(testarray)
print(a, ae * np.sqrt(2 * pt580[2]), b, be * np.sqrt(2 * pt580[2]), c, ce * np.sqrt(2 * pt580[2]), d, de * np.sqrt(2 * pt580[2]), e, ee * np.sqrt(2 * pt580[2]), chi2 / (2 * pt580[2]))
testarrayre = np.real(testarray)
testarrayim = np.imag(testarray)
resre = []
resres = []
resim = []
resims = []
for j in range(24):
    _, _, t = AutoCorrelationSingleVariable(testarrayre[:, j], s=2.0)
    v, s = JacknifeMean(testarrayre[:, j])
    resre.append(v)
    resres.append(s * np.sqrt(2 * t))
    _, _, t = AutoCorrelationSingleVariable(testarrayim[:, j], s=2.0)
    v, s = JacknifeMean(testarrayim[:, j])
    resim.append(v)
    resims.append(s * np.sqrt(2 * t))
print(PrintAsMatlabArray(resre))
print(PrintAsMatlabArray(resres))
print(PrintAsMatlabArray(resim))
print(PrintAsMatlabArray(resims))
plt.plot(resre, resim)
xx = np.array([m * 0.1 for m in range(240)])
yy = polyakov_fit580b2(xx, [a, b, c, d, e])
yyre = np.real(yy)
yyim = np.imag(yy)
plt.plot(yyre, yyim)
plt.show()


"""
0.2872060294245133 0.004075796077113506 -0.00672455072566631 0.002806639984584521 0.004089039260148491 0.0024928468889779246 -0.023609489095259898 0.0063854188718763015 0.01794911880344403 0.0049638437658799765 0.03858969148469004
[0.28123737,0.27492571,0.26863995,0.25964397,0.25510025,0.25566402,0.26429109,0.28073872,0.29943841,0.31456431,0.32616166,0.33465192,0.33786092,0.33379387,0.32196537,0.30717643,0.28664307,0.27361382,0.26280532,0.25752794,0.26183652,0.26951622,0.27668357,0.28208351]
[0.0092813,0.00953128,0.01082616,0.01033758,0.01294405,0.01135748,0.00864858,0.00517612,0.00338796,0.00351672,0.00294202,0.0032058,0.00317495,0.00370388,0.00336968,0.00426069,0.00467689,0.00459145,0.00431783,0.00347377,0.00381417,0.00501713,0.00661501,0.00784187]
[1.49052224e-03,3.86713933e-03,8.57057780e-03,1.07192044e-02,1.26199548e-02,1.14037732e-02,9.22036453e-03,5.53996441e-03,2.17787125e-03,8.69771773e-05,2.89364699e-04,1.22375354e-03,9.07578848e-04,-1.06215285e-04,6.99303597e-04,-2.58230302e-04,-1.93391907e-03,-3.35969199e-03,-5.62233887e-03,-8.30924914e-03,-7.30741062e-03,-7.32970728e-03,-5.83938452e-03,-2.65171583e-03]
[0.00112884,0.00108369,0.00120277,0.00141127,0.00161603,0.00170298,0.00117956,0.00140493,0.00112831,0.00094661,0.00088841,0.00085452,0.00090305,0.00093298,0.00102351,0.00109323,0.00104009,0.00151839,0.00143105,0.00133042,0.00151838,0.00120679,0.00120813,0.00120237]
0.29762606134953357 0.003049473463629256 0.0030380467325575038 0.0020428791801190235 0.002074678036388598 0.0018880881420864216 0.01441197042887903 0.003976742657250897 0.004043705131199402 0.00308649208706474 0.1159627085522742
[0.31479357,0.30889054,0.3008019,0.29077368,0.28884048,0.28767257,0.29031083,0.29080927,0.29074756,0.29568553,0.30814192,0.32400432,0.32637979,0.31711164,0.30314982,0.29245302,0.28535296,0.28097465,0.27823057,0.27496938,0.2734318,0.27942775,0.29531128,0.31013968]
[0.00585117,0.00521604,0.00539384,0.00510343,0.00420707,0.00442639,0.00384513,0.0027164,0.0037226,0.0031728,0.00361034,0.00323534,0.00229177,0.00239638,0.0024993,0.00318346,0.00357305,0.00435387,0.00532663,0.00749555,0.008779,0.00964302,0.00885333,0.00681291]
[-3.02138499e-05,-8.61075522e-04,-2.47891700e-03,-6.42760756e-03,-5.58828313e-03,-5.09362183e-03,-1.02280628e-03,3.91908723e-03,5.25920544e-03,2.44773902e-03,1.37834556e-03,-1.01599997e-03,8.57522014e-04,6.33909372e-04,-5.36709696e-04,-1.45436160e-03,-1.97974614e-03,-3.81820227e-04,3.15433174e-03,5.53336565e-03,7.07177960e-03,2.97305066e-03,1.65078551e-03,-2.90504021e-04]
[0.00098906,0.00088389,0.00105529,0.001401,0.00107572,0.00103514,0.00112749,0.00115136,0.00103783,0.0010614,0.00098096,0.00091614,0.0009122,0.00092639,0.00099382,0.00093909,0.00114099,0.00114539,0.00140696,0.00129368,0.00143369,0.00117613,0.001088,0.00097332]
"""
