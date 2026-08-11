from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads = "H:\\BS01\\BS01Nt6\\Chiral\\"


lst = ["530", "535", "540", "545", "550", "555", "560", "565", "570", "575", "580"]

cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(len(lst)):
    fileNames = folderHeads + "BS01__{}_condensatepCCLightChiralKS.csv".format(lst[i])
    testarray = LoadMathematicaCSV(fileNames)
    testarrayre = np.real(testarray)
    v, s = JacknifeMean(testarrayre, lst[i])
    _, _, t = AutoCorrelationSingleVariable(testarrayre)
    cv.append(-v)
    sv.append(s * np.sqrt(2 * t))
    tvalues.append(t)
    v, s = JacknifeCumulant(testarrayre, lst[i])
    print("==============", i, "with {} results".format(len(testarrayre)))
    print(v)
    print(s)
    values.append(v)
    stds.append(s * np.sqrt(2 * t))

print(PrintAsMathematicaArray(cv, "v"))
print(PrintAsMathematicaArray(sv, "s"))
print(PrintAsMathematicaArray(tvalues, "t"))
print(PrintAsMathematicaArray(values, "susp"))
print(PrintAsMathematicaArray(stds, "suspe"))


errorbar(range(len(values)), [cv], [sv])
errorbar(range(len(values)), [values], [stds])

"""
v={0.42694061,0.37479436,0.31777675,0.26370413,0.21400646,0.17388158,0.15536841,0.14388803,0.13537655,0.12889493,0.12351807};
s={3.27888695*^-04,5.10113729*^-04,5.94142878*^-04,4.51905701*^-04,2.37322696*^-03,1.23027545*^-03,4.51374584*^-04,1.20218280*^-04,9.53289115*^-05,8.58464632*^-05,4.22586433*^-05};
t={10.85199668,19.67790558,24.13587039,18.57054297,146.8973696,102.54862005,42.94604317,9.18567727,8.15294456,9.09225058,3.38730263};
susp={6.93491713*^-06,9.25663615*^-06,1.02380410*^-05,7.69784352*^-06,2.68387675*^-05,1.03317272*^-05,3.32084867*^-06,1.10135639*^-06,7.80248277*^-07,5.67376656*^-07,3.69041439*^-07};
suspe={1.18358127*^-06,2.17123994*^-06,3.00423126*^-06,1.59456910*^-06,1.45128046*^-05,5.67483828*^-06,1.11819002*^-06,1.76975076*^-07,1.26681563*^-07,8.84931567*^-08,3.65669028*^-08};

susp 量纲: 1111.1710
chi = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 Nx^3 Nt)
chi / T^2 = 0

可以看到，这里要求chi ~ -2，D~1，然而D~1是不可能的。我们考虑3种定义：

L = \int dx psibar D psi 积分 这里D~m~-1
L = psibar D psi  矩阵 这里D~3，因为L~0

所以，实际上的定义应为：
chi = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/V4 Ntaste^2 =  Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 a^4 Nx^3 Nt)
chi ~ -2，是因为D~3
我们以D=a^4 Dhat, d = 2a Dhat来表示，其中d是实际测量的。那么d=2 a^-3 D  2 a^-3*d^-1 = D^-1

(<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2是实际测量的
4 * a^-6 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2 = Nx^3Nt * (<tr[D-1]^2> - <tr[D-1]>^2)/(Nx^3 Nt)^2 = (<tr[D-1]^2> - <tr[D-1]>^2)/(Nx^3 Nt)
4 * a^4 * a^-6 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2 = (<tr[D-1]^2> - <tr[D-1]>^2)/(a^4 Nx^3 Nt)
4/16 * 4 * a^4 * a^-6 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2 = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 a^4 Nx^3 Nt)
左边整理一下
a^-2 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2 = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 a^4 Nx^3 Nt)

a^-2 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2/(a^-2/Nt^2) = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 a^4 Nx^3 Nt)/T^2
Nt^2 * Nx^3Nt * (<tr[d-1]^2> - <tr[d-1]>^2)/(Nx^3Nt)^2 = Nf^2 (<tr[D-1]^2> - <tr[D-1]>^2)/(16 a^4 Nx^3 Nt)/T^2


"""