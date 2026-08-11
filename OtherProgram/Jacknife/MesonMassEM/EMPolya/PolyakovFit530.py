from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

from MesonMassEM.Mass.CorrelatedFitGPT import PolyaFit, polyakov_fit, PolyaFit2, polyakov_fit2

hd = "530"
folderHeads = f"H:\\MesonMassElectric\\{hd}01\\Polyakov\\"
# folderHeads = "H:\\BS01Polyakov\\"

pt530 = [0.6029984304961401, 0.5965243156398604, 0.5784913864290465, 0.65875253976466, 0.589699822873557, 0.5671786190503038, 0.5913288100050554, 0.6535853729622049]

cv = []
sv = []
values = []
stds = []
tvalues = []


fileNames = folderHeads + f"EC{hd}_01_M0__1_polyakov_XSlice.csv"
testarray = LoadMathematicaCSV(fileNames)[100:, :]
a, ae, b, be, c, ce, chi2 = PolyaFit(testarray)
print(a, ae * np.sqrt(2 * pt530[1]), b, be * np.sqrt(2 * pt530[1]), c, ce * np.sqrt(2 * pt530[1]), chi2 / (2 * pt530[1]))
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
yy = polyakov_fit(xx, [a, b, c])
yyre = np.real(yy)
yyim = np.imag(yy)
plt.plot(yyre, yyim)
plt.show()

fileNames = folderHeads + f"EC{hd}_01_M0__2_polyakov_XSlice.csv"
testarray = LoadMathematicaCSV(fileNames)
a, ae, b, be, c, ce, chi2 = PolyaFit2(testarray)
print(a, ae * np.sqrt(2 * pt530[2]), b, be * np.sqrt(2 * pt530[2]), c, ce * np.sqrt(2 * pt530[2]), chi2 / (2 * pt530[2]))
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
yy = polyakov_fit2(xx, [a, b, c])
yyre = np.real(yy)
yyim = np.imag(yy)
plt.plot(yyre, yyim)
plt.show()


"""
-3.2188612502585764e-05 0.00020024615394190532 -0.002931447840679531 0.00020490668743692046 0.0027408435333646118 0.0002026636057951175 0.5930708990362105
[-0.0002574,-0.00109914,-0.00247389,-0.00309849,-0.00270931,-0.00352699,-0.00269293,-0.00175113,-0.00113063,0.00229136,0.00361799,0.00563479,0.00623929,0.00284728,0.00568721,0.00195158,0.00168726,-0.00044804,-0.00233455,-0.00305005,-0.00354581,-0.00253175,-0.0005495,0.00027239]
[0.00102722,0.00096533,0.00100523,0.00097018,0.00108271,0.00093876,0.00101351,0.00104135,0.00103521,0.00110781,0.00121189,0.00104579,0.0011241,0.00104783,0.00098761,0.0009692,0.00115509,0.00106829,0.00097235,0.00101926,0.00091956,0.00096503,0.00103769,0.00100781]
[1.63234411e-03,1.55142904e-03,4.38459898e-03,3.93905247e-03,6.58102173e-03,4.87171841e-03,3.37024112e-03,3.65535136e-04,-9.25882702e-04,-2.54630297e-03,-5.49830577e-04,1.41427921e-04,1.03157535e-03,1.57023517e-05,3.03770418e-03,2.89954777e-05,7.34419879e-06,-2.88858867e-03,-2.32552090e-03,-3.14519251e-03,-4.25176623e-03,-5.06816107e-03,-3.37634314e-03,-4.96393687e-04]
[0.00103534,0.00099585,0.00096512,0.00099171,0.00099948,0.00103044,0.00107682,0.00120073,0.00107499,0.00100296,0.00093537,0.00099645,0.00097064,0.00097251,0.00106114,0.00092609,0.00099815,0.00105176,0.00125036,0.00098943,0.00109784,0.00096458,0.00110585,0.00111051]
7.329547017498246e-06 0.0001868227491690899 0.002958547396978565 0.00019133869275884016 0.0030863761289669895 0.0001809559209385577 0.5153650566428165
[0.00665205,0.00375947,0.00149035,-0.00249099,-0.00285938,-0.0013229,0.00031783,-0.00118266,-0.00258043,-0.00472579,-0.00100488,0.00390665,0.00563337,0.00316174,-0.00018291,-0.00490367,-0.00197415,-0.00123958,-0.00064375,-0.00096902,-0.00256307,-0.00224677,0.00200797,0.00485783]
[0.00104458,0.00103972,0.00107967,0.00094534,0.00097297,0.00104165,0.00105484,0.00111618,0.00104171,0.00093709,0.00094189,0.00122319,0.00107597,0.00098356,0.00102356,0.00104244,0.00101398,0.00102689,0.00099874,0.00094808,0.00105789,0.00101524,0.00092937,0.00101139]
[0.0004626,0.00076925,0.00056202,-0.00211843,-0.00503598,-0.00303664,-0.00039585,0.00453084,0.00585085,0.00396312,0.00017911,0.00110087,0.0013547,0.00141629,0.00029806,-0.00369219,-0.00614605,-0.00369773,-0.00022359,0.00371659,0.00535903,0.00219137,-0.00107345,-0.00325656]
[0.00095043,0.0010199,0.00093982,0.00111747,0.00108021,0.00104173,0.00095531,0.00095874,0.00097379,0.00100881,0.00105637,0.00098731,0.00104362,0.00090706,0.0010176,0.00101102,0.00097417,0.0010878,0.00101213,0.00101203,0.00111948,0.00104013,0.00094934,0.00095935]

"""
