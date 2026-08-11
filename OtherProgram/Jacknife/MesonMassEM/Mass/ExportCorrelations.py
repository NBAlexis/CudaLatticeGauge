"""
我们检查：

du=ud?

2=3?

9=10?

输出：

f_id_electric

f = 0,1,2 代表(dd, uu, du)
id = 0-7 代表 [[0], [11], [2, 3], [4], [8], [9, 10]]
electric = 0-7，代表外加电场强度

"""
import numpy as np
from matplotlib import pyplot as plt

from JacknifePrograms import LoadMathematicaCSV

header = "53001"
folderHeads = "H:\\MesonMassElectric\\" + header + "\\Meson\\"
flst = [["uu"], ["dd"], ["du", "ud"]]
typelst = [[0], [11], [2, 3], [4], [8], [9, 10]]

def csymtry(c):
    sizeofc = np.shape(c)
    retv = np.zeros((sizeofc[0], sizeofc[1] // 2 + 1))
    for i in range(1, sizeofc[1] // 2):
        retv[:, i] = 0.5 * (c[:, i] + c[:, sizeofc[1] - i])
    retv[:, 0] = c[:, 0]
    retv[:, sizeofc[1] // 2] = c[:, sizeofc[1] // 2]
    return retv

for k in range(8):
    for i in range(len(flst)):
        for j in range(len(typelst)):
            correlator_arrays = []
            for ii in range(len(flst[i])):
                for jj in range(len(typelst[j])):
                    fileNames = folderHeads + f"GFEC530_01_M0_X_{flst[i][ii]}{typelst[j][jj]}_{k}.csv"
                    testarray = np.real(LoadMathematicaCSV(fileNames))
                    correlator_arrays.append(testarray)
            """
            当flst[i], typelst[j]有多个的时候，把它们画到一起
            """
            res = correlator_arrays[0]
            if len(correlator_arrays) > 1:
                """
                for ii in range(len(correlator_arrays)):
                    v = np.mean(correlator_arrays[ii], axis=0)
                    s = np.var(correlator_arrays[ii], axis=0)
                    plt.errorbar(range(len(v)), v, yerr=s, fmt="o")
                plt.title(f"i={i},j={j},k={k}")
                plt.show()
                """
                for ii in range(1, len(correlator_arrays)):
                    res = res + correlator_arrays[ii]
            res = csymtry(res)
            np.savetxt(f'../Data/{header}_f{i}_t{j}_e{k}.csv', res, delimiter=',')
            print(f'../Data/{header}_f{i}_t{j}_e{k}.csv saved')
