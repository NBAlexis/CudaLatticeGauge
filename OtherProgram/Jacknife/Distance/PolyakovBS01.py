from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads = "H:\\BS01\\BS01Nt6\\Polyakov\\"


lst = ["530", "535", "540", "545", "550", "555", "560", "565", "570", "575", "580"]

cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(len(lst)):
    fileNames = folderHeads + "BS01__{}_polyakov.csv".format(lst[i])
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
v={0.00341303,0.0038655,0.00584383,0.01047444,0.03277648,0.06641688,0.08267707,0.09421701,0.10320332,0.11080497,0.11884657};
s={5.99201737*^-05,7.84350908*^-05,1.77536052*^-04,4.38926008*^-04,3.72616999*^-03,1.96046182*^-03,9.71153141*^-04,3.24071873*^-04,2.89746503*^-04,6.23954976*^-04,3.89787527*^-04};
t={0.8714636,1.15827743,3.66294882,13.10916217,151.15646985,84.46133476,31.1694381,5.67451739,4.69285755,14.7620313,8.19474625};
susp={2.88399775*^-06,3.71797320*^-06,6.02338058*^-06,1.02874026*^-05,6.42978759*^-05,3.18534794*^-05,2.11809047*^-05,1.29554286*^-05,1.25226740*^-05,1.84611361*^-05,1.29783178*^-05};
suspe={1.43157622*^-07,2.07435424*^-07,5.45819463*^-07,1.80855573*^-06,3.69722735*^-05,1.70855443*^-05,6.20404107*^-06,1.54718597*^-06,1.47472786*^-06,3.59060667*^-06,2.07036972*^-06};

susp的量纲

测量的是无量纲的P=\int d^3x p / V3
需要的是 chi=Nx^3(<P>^2-<P^2>)

但chi/T^3 = 0，所以，这里P应该是\int d^3 x p ~3，


"""