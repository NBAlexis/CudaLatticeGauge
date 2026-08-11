import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *
import matplotlib.pyplot as plt

from MesonMassEM.Mass.UsefulFunctions import jackknife_meff, staggered_correlated_fit, jackknife_meff2

"""
dd-0 means:

0 means: x,y,z,t
dd means dirac operator
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

cv = []
sv = []
values = []
stds = []
tvalues = []

t0_0=[1.8951564773357852, 1.8077776285910017, 2.1946299065815835, 1.990912303869544, 2.079855920214655, 2.5339498646686325, 1.9359976596659145, 1.7001348019984017]
t0_1=[0.5, 0.5006775446534039, 0.5, 0.5035645932041654, 0.500410055055439, 0.5015858973486312, 0.5022722705371013, 0.5005470051010357]
t0_2=[0.5, 0.5093389034885256, 0.501093118503719, 0.5, 0.5007772148515922, 0.5, 0.5021374216384019, 0.5]
t0_3=[0.5, 0.5057604991624196, 0.5, 0.5, 0.5, 0.5, 0.500601505613508, 0.5001208717234339]
t0_4=[0.5025920683666747, 0.5042377795198644, 0.5, 0.5008509246634495, 0.5037768147535879, 0.5014176646269032, 0.5008449891192723, 0.5006369704132312]
t0_5=[0.5010103662061772, 0.5, 0.5048361330271944, 0.5029681942065802, 0.5, 0.5009024029430554, 0.5, 0.5015919597170022]
t1_0=[1.8951564773357852, 1.9855314099039842, 2.2630674063074196, 2.017871206484635, 1.828738872347699, 2.4314211165244823, 2.2580996646621387, 2.4495390740346696]
t1_1=[0.5, 0.5, 0.5007967863609972, 0.5, 0.5011841819487235, 0.5, 0.5008447252021568, 0.5]
t1_2=[0.5, 0.5011265383606144, 0.5, 0.5057217360292533, 0.5, 0.5009088513127982, 0.5001126201018703, 0.5295632409638937]
t1_3=[0.5, 0.5017670165715111, 0.5, 0.5, 0.5068143906748243, 0.5, 0.5004754991140478, 0.5025053475580304]
t1_4=[0.5025920683666747, 0.5, 0.5051326311884138, 0.5007175748200305, 0.5010454904988729, 0.5, 0.5028526658250584, 0.5]
t1_5=[0.5010103662061772, 0.5, 0.5009946223881415, 0.5, 0.5002983724490149, 0.5000261658127404, 0.5014506098493993, 0.5317288463889387]
t2_0=[1.8951564773357852, 1.895264902007808, 2.4833655233535312, 2.117904731451813, 1.976755391480002, 2.6287616011817017, 2.135726720459151, 2.3080403836883567]
t2_1=[0.5, 0.5, 0.5005331227293764, 0.5002897001369178, 0.5, 0.5009558867080589, 0.5, 0.5]
t2_2=[0.5, 0.5004507662443937, 0.5016887711255547, 0.505394134553677, 0.501503519123976, 0.5054397248320796, 0.5012184193731878, 0.5]
t2_3=[0.5, 0.5009185638895226, 0.5, 0.5008995340733087, 0.5, 0.5, 0.5032329401389063, 0.5]
t2_4=[0.5025920683666747, 0.5158195190931527, 0.5002146738542255, 0.5202130471176772, 0.5, 0.5, 0.5, 0.5]
t2_5=[0.5010103662061729, 0.5, 0.5018538754887764, 0.5, 0.5, 0.5005057400906975, 0.5021138338126407, 0.5029163897011946]

pmlst = [1, 2, 1, 1, 2, 2]

tp0t = [t0_0, t0_1, t0_2, t0_3, t0_4, t0_5]
tp1t = [t1_0, t1_1, t1_2, t1_3, t1_4, t1_5]
tp2t = [t2_0, t2_1, t2_2, t2_3, t2_4, t2_5]
tfpt = [tp0t, tp1t, tp2t]

all_eff_mass = []

for f in range(3):
    for tp in range(6):
        for i in range(8):
            c = np.loadtxt(f"../Data/{hd}_f{f}_t{tp}_e{i}.csv", delimiter=",")
            if 0 == i:
                c = c[100:, :]
            if 1 == pmlst[tp]:
                c = cp(c)
            elif 2 == pmlst[tp]:
                c = cm(c)
            # print(c)
            tlst = tfpt[f][tp]
            if tp >= 3:
                vp, sp, xa = jackknife_meff2(c)
                all_eff_mass.append(vp)
                all_eff_mass.append(sp * np.sqrt(2 * tlst[i]))
                plt.errorbar(xa, vp, yerr=sp)
            else:
                vp, sp, xa = jackknife_meff(c)
                all_eff_mass.append(vp)
                all_eff_mass.append(sp * np.sqrt(2 * tlst[i]))
                plt.errorbar(xa, vp, yerr=sp)
        if tp >= 3:
            plt.ylim([0.3, 2])
        plt.show()

all_eff_mass = np.array(all_eff_mass)
np.savetxt(f"../DrawData/{hd}_effmass.csv", all_eff_mass, delimiter=",")

"""
for i in range(8):
    fileNames = folderHeads + f"GFEC530_01_M0_X_uu11_{i}.csv"
    testarray = np.real(LoadMathematicaCSV(fileNames))
    plt.plot(np.mean(testarray, axis=0))
    plt.show()
    csymmetric = csymtry(testarray)
    print(np.mean(csymmetric, axis=0))
    cp_thischannel = cp(csymmetric)
    cm_thischannel = cm(csymmetric)
    print(np.mean(cp_thischannel, axis=0))
    print(np.mean(cm_thischannel, axis=0))
    vp, sp, xa = jackknife_meff(cp_thischannel)
    vm, sm, _ = jackknife_meff(cm_thischannel)
    plt.errorbar(xa, vp, yerr=sp)
    plt.errorbar(xa, vm, yerr=sm)
    plt.show()

    # res = staggered_correlated_fit(csymmetric, 24, 2, False)
    # print(f"m+={res["params"][1]}+-{res["params_err"][1]}, m-={res["params"][3]}+-{res["params_err"][3]}, chi2/dof={res["chi2_dof"]}")
    mplus_mean, Aplus_mean, mminus_mean, Aminus_mean, dmplus, dAplus, dmminus, dAminus, chi2_dof = staggered_correlated_fit(csymmetric, 24, 5, 2)
    print(f"m+={mplus_mean}+-{dmplus}, m-={mminus_mean}+-{dmminus}, chi2/dof={chi2_dof}")
    # print(np.mean(testarray, axis=0))
"""

# print(PrintAsMathematicaArray(values, "susp"))
# print(PrintAsMathematicaArray(stds, "suspe"))
# print(PrintAsMathematicaArray(tvalues, "t"))

# errorbar(range(len(values)), [cv], [sv])
# errorbar(range(len(values)), [values], [stds])

