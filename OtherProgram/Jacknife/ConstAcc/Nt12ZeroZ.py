import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from ConstAcc.UsefullFunction import MeasureValueType
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean
from Visualization import errorbar

headnames = ["68", "69", "70", "76", "80"]
glst = ["00", "03", "06", "09", "12", "15", "18", "21", "24", "27", "30"]

mesureTypeList = [
    MeasureValueType.POLYA,
    MeasureValueType.CHIRAL,
    MeasureValueType.G3,
    MeasureValueType.G4
]
filepatternList = [
"H:\\ConstAccLarge\\ACCQM{}\\Polyakov\\ACCQM{}__{}_polyakov_ZSlice.csv",
"H:\\ConstAccLarge\\ACCQM{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightChiralKS.csv",
"H:\\ConstAccLarge\\ACCQM{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma3.csv",
"H:\\ConstAccLarge\\ACCQM{}\\Chiral\\ACCQM{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv"
]

ylabellst = ["$P$", "$c$", "$c_3$", "$Im(c_4)$"]
savenames = ["polya", "chiral", "gamma3", "gamma4"]

tomeasure = 3
measureType = mesureTypeList[tomeasure]
filepattern = filepatternList[tomeasure]

todraw = []
todrawe = []

for hid in range(len(headnames)):
    vlst = []
    elst = []
    tlst = []
    for gid in range(len(glst)):
        filetoload = filepattern.format(headnames[hid], headnames[hid], glst[gid])
        arr = LoadMathematicaCSV(filetoload)
        arr = arr[:, 3]
        if MeasureValueType.POLYA == measureType:
            arr = np.abs(arr) / 3
        elif MeasureValueType.CHIRAL == measureType:
            """
            arr was -Tr[D]
            we need Tr[D] which is -real(arr)
            """
            arr = -np.real(arr) / 2
        elif MeasureValueType.G3 == measureType:
            """
            arr was -Tr[D]
            we need i * (-Tr[D]) which is -i Tr[D] which is -imag(arr)
            """
            arr = -np.imag(arr) / 2
        elif MeasureValueType.G4 == measureType:
            """
            arr was -Tr[D]
            we need Im[-Tr[D]] which is i * Tr[D] which is imag(arr)
            """
            arr = np.imag(arr) / 2
        v, s = JacknifeMean(arr)
        _, _, t = AutoCorrelationSingleVariable(arr)
        vlst.append(v)
        elst.append(s)
        tlst.append(t)
    vlst = np.array(vlst)
    elst = np.array(elst)
    tlst = np.array(tlst)
    todraw.append(vlst)
    todrawe.append(elst * np.sqrt(2 * tlst))

todraw = np.array(todraw)
todrawe = np.array(todrawe)
ledlst = ["T=278.1 MeV", "T=289.1 MeV", "T=303.0 MeV", "T=349.8 MeV", "T=383.6 MeV"]
xlst = [0.03 * n * 12 / 3.1415926 for n in range(11)]

errorbar(xlst, todraw, todrawe, xlabel='$g/(T\pi)$', ylabel=ylabellst[tomeasure], legends=ledlst, savefile="fig" + savenames[tomeasure] + ".eps")
