import numpy as np
from matplotlib import pyplot as plt

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant

# foldername = "I:/BetaGradient/BetaScan/BS529536/Polyakov/"
foldername = "H:/MesonMassElectric/Chrial/"
foldername2 = "H:/MesonMassElectric/Polyakov/"
betas = ["490", "550"]
fermiontype = ["Heavy", "Light"]
condensationtype = ["ChiralKS",
                    "CMTKSGamma4", "CMTKSGamma5",
                    "CMTKSSigma12", "CMTKSSigma34"]

factors = [-1, 1j, -1, -1, -1]
maxz = 12
"""
490:
re chiral
re sigma12
re sigma34

550:
re gamma5

"""
drawchiral = False
drawpolya = True

for b in betas:
    for c in range(len(condensationtype)):
        hrelst = []
        lrelst = []
        htlst = []
        ltlst = []
        hslst = []
        lslst = []
        hcumlst = []
        hcumslst = []
        lcumlst = []
        lcumslst = []
        for f in range(len(fermiontype)):
            for em in range(11):
                fileNames = foldername + "{}01\\EC{}\\EC{}_01__{}_condensatepCC{}{}.csv".format(b, em, b, em, fermiontype[f], condensationtype[c])
                testarray = LoadMathematicaCSV(fileNames)
                arr = np.real(factors[c] * testarray)
                v, s = JacknifeMean(arr)
                _, _, t = AutoCorrelationSingleVariable(arr)
                s = s * np.sqrt(2 * t)
                if 0 == f:
                    hrelst.append(v)
                    htlst.append(t)
                    hslst.append(s)
                else:
                    lrelst.append(v)
                    ltlst.append(t)
                    lslst.append(s)
                if 0 == c:
                    v, s = JacknifeCumulant(arr)
                    s = s * np.sqrt(2 * t)
                    if 0 == f:
                        hcumlst.append(v)
                        hcumslst.append(s)
                    else:
                        lcumlst.append(v)
                        lcumslst.append(s)
        if drawchiral:
            plt.errorbar(range(11), hrelst, yerr=hslst)
            plt.errorbar(range(11), lrelst, yerr=lslst)
            plt.title("{}-{}".format(b, condensationtype[c]))
            plt.show()
            if 0 == c:
                plt.errorbar(range(11), hcumlst, yerr=hcumslst)
                plt.errorbar(range(11), lcumlst, yerr=lcumslst)
                plt.title("{}-{}cum".format(b, condensationtype[c]))
                plt.show()
        if 0 == c:
            np.savetxt(foldername + "{}-{}.csv".format(b, condensationtype[c]), np.vstack((
                np.array([hrelst]),
                np.array([hslst]),
                np.array([hcumlst]),
                np.array([hcumslst]),
                np.array([htlst]),
                np.array([lrelst]),
                np.array([lslst]),
                np.array([lcumlst]),
                np.array([lcumslst]),
                np.array([ltlst]))),
                delimiter=',')
        else:
            np.savetxt(foldername + "{}-{}.csv".format(b, condensationtype[c]), np.vstack((
                np.array([hrelst]),
                np.array([hslst]),
                np.array([htlst]),
                np.array([lrelst]),
                np.array([lslst]),
                np.array([ltlst]))),
                delimiter=',')
    abslst = []
    slst = []
    tlst = []
    cumlst = []
    cumslst = []
    zdata = []
    for em in range(11):
        fileNames = foldername2 + "{}01\\EC{}_01__{}_polyakov.csv".format(b, b, em)
        testarray = LoadMathematicaCSV(fileNames)
        arr = np.abs(testarray)
        v, s = JacknifeMean(arr)
        _, _, t = AutoCorrelationSingleVariable(arr)
        s = s * np.sqrt(2 * t)
        abslst.append(v)
        slst.append(s)
        tlst.append(t)
        v, s = JacknifeCumulant(arr)
        s = s * np.sqrt(2 * t)
        cumlst.append(v)
        cumslst.append(s)
        absz = []
        rez = []
        imz = []
        rezs = []
        imzs = []
        abszs = []
        fileNames = foldername2 + "{}01\\EC{}_01__{}_polyakov_ZSlice.csv".format(b, b, em)
        arr = LoadMathematicaCSV(fileNames)
        for z in range(0, maxz):
            thisZSlice = arr[:, z]
            thisZSlice = thisZSlice.flatten()
            aarr = np.abs(thisZSlice)
            rarr = np.real(thisZSlice)
            iarr = np.imag(thisZSlice)
            v, s = JacknifeMean(aarr)
            absz.append(v)
            abszs.append(s * np.sqrt(2 * t))
            v, s = JacknifeMean(rarr)
            rez.append(v)
            rezs.append(s * np.sqrt(2 * t))
            v, s = JacknifeMean(iarr)
            imz.append(v)
            imzs.append(s * np.sqrt(2 * t))
        if drawpolya:
            # plt.errorbar(range(12), absz, yerr=abszs)
            plt.errorbar(range(12), rez, yerr=rezs)
            plt.errorbar(range(12), imz, yerr=imzs)
            plt.title("{}-{}".format(b, em))
            plt.show()
        zdata.append(absz)
        zdata.append(abszs)
        zdata.append(rez)
        zdata.append(rezs)
        zdata.append(imz)
        zdata.append(imzs)
    np.savetxt(foldername2 + "p{}.csv".format(b), np.vstack((
        np.array([abslst]),
        np.array([slst]),
        np.array([cumlst]),
        np.array([cumslst]),
        np.array([tlst])
    )), delimiter=',')
    np.savetxt(foldername2 + "p{}-z.csv".format(b), np.array(zdata), delimiter=',')


