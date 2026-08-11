import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

from MesonMassEM.Mass.UsefulFunctions import jackknife_meff, staggered_correlated_fit, jackknife_meff2

"""
we calculate t_int
"""

hd = "53001"

def cp(c):
    sizeofc = np.shape(c)
    retv = np.zeros((sizeofc[0], sizeofc[1] - 1))
    for i in range(sizeofc[1] - 1):
        retv[:, i] = c[:, i] + c[:, i + 1]
    return retv

def cm(c):
    sizeofc = np.shape(c)
    retv = np.zeros((sizeofc[0], sizeofc[1] - 1))
    for i in range(sizeofc[1] - 1):
        retv[:, i] = (c[:, i] - c[:, i + 1]) * ((-1)**i)
    return retv

def meffv(c):
    sizeofc = np.shape(c)
    retv = np.zeros((sizeofc[0], sizeofc[1] - 2))
    for i in range(sizeofc[1] - 2):
        retv[:, i] = (c[:, i] + c[:, i + 2]) / (2 * c[:, i + 1])
    return retv

cv = []
sv = []
values = []
stds = []
tvalues = []

pmlst = [1, 2, 1, 1, 2, 2]

all_eff_mass = []
tallcase = []
for f in range(3):
    for tp in range(6):
        tlst = []
        for i in range(8):
            c = np.loadtxt(f"../Data/{hd}_f{f}_t{tp}_e{i}.csv", delimiter=",")
            if i != 0:
                c = c[100:, :]
            if 1 == pmlst[tp]:
                c = cp(c)
            elif 2 == pmlst[tp]:
                c = cm(c)
            # c = cp(c)
            # print(c)
            mef = meffv(c)
            mef = np.mean(mef, axis=1)
            # mef = np.acosh(mef)
            # print(np.shape(mef))
            # print(mef)
            # v, s, t = AutoCorrelationSingleVariable(mef[:, 0])
            v, s, t = AutoCorrelationSingleVariable(mef, s=2.0)
            tlst.append(float(t))
        print(f"t{f}_{tp}={tlst}")
        # if tp >= 3:
        #     plt.ylim([0.3, 2])
        # plt.show()

# all_eff_mass = np.array(all_eff_mass)
# np.savetxt(f"../DrawData/{hd}_effmass.csv", all_eff_mass, delimiter=",")



# print(PrintAsMathematicaArray(values, "susp"))
# print(PrintAsMathematicaArray(stds, "suspe"))
# print(PrintAsMathematicaArray(tvalues, "t"))

# errorbar(range(len(values)), [cv], [sv])
# errorbar(range(len(values)), [values], [stds])

