from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

# folderHeads = "G:\\BS01Nt12Chiral\\"
folderHeads = "F:\\Builds\\FloatBuild-04-30-sm86\\Release\\"


# lst = ["530", "535", "540", "540", "545", "550", "555", "560", "565", "570", "575", "580"]
# lst = ["63", "64", "65", "66", "67", "68", "69", "70", "71", "72"]
# lst = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
lst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]
# lst2 = ["BSQ", "BSQ2", "BSQ3", "BSQ4", "BSQ5", "BSQ6", "BSQ7"]
# lst = ["51400", "51775", "51900", "52150", "52400", "52525", "52900", "53275", "53400", "53650", "53900", "54025", "54400"]
# lst = ["52100", "52475", "52600", "52850", "53100", "53225", "53600", "53975", "54100", "54350", "54600", "54725", "55100"]
# lst2 = ["0-30", "30-60"]
lst2 = [""]
# lst2 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "T", "S", "T", "U"]

cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(len(lst)):
    for j in range(len(lst2)):
        fileNames = folderHeads + "BSQ__{}_polyakov.csv".format(lst[i])
        # fileNames = folderHeads + "BS01__{}_condensatepCCLightChiralKS.csv".format(lst[i])

        testarray = LoadMathematicaCSV(fileNames)
        if 0 == j:
            testarrayre = np.abs(testarray)
            # testarrayre = np.real(testarray)
        else:
            # testarrayre = np.hstack((testarrayre, np.abs(testarray)))
            testarrayre = np.hstack((testarrayre, np.real(testarray)))
    # testarrayre = testarrayre * testarrayre
    v, s = JacknifeMean(testarrayre, lst[i])
    _, _, t = AutoCorrelationSingleVariable(testarrayre)
    cv.append(-v)
    sv.append(s * np.sqrt(2 * t))
    v, s = JacknifeCumulant(testarrayre, lst[i])
    print("==============", i, "with {} results".format(len(testarrayre)))
    print(v)
    print(s)
    # print(t)
    values.append(v)
    stds.append(s * np.sqrt(2 * t))
    tvalues.append(t)


print(PrintAsMathematicaArray(cv, "chiral"))
print(PrintAsMathematicaArray(sv, "chirale"))
print(PrintAsMathematicaArray(values, "susp"))
print(PrintAsMathematicaArray(stds, "suspe"))
print(PrintAsMathematicaArray(tvalues, "t"))

errorbar(range(len(values)), [cv], [sv])
errorbar(range(len(values)), [values], [stds])

"""
lst = ["52", "53", "54", "55", "56", "57", "58"]
16*16*16*6 101-1000, am=0.1
chiral={0.5156876,0.44909172,0.35328924,0.2528227,0.17674085,0.13983815,0.12605397};
chirale={0.00042646,0.0009435,0.00087268,0.00093335,0.003174,0.0004118,0.00017358};
susp={1.70144451*^-05,2.33300490*^-05,2.57166656*^-05,2.30320767*^-05,4.94008313*^-05,5.40654808*^-06,2.09851233*^-06};
suspe={2.45011017*^-06,5.74045064*^-06,5.92036481*^-06,6.53369178*^-06,2.94187595*^-05,1.20747389*^-06,3.66540188*^-07};


lst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]
24 24 24 6, 101-3000, quenching. Polyakov


"""