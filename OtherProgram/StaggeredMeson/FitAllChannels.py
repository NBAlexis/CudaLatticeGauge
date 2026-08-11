import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LinearLocator

from UsefulFunctions import fitcorrelation

nt = 64
t_start = 10
# spacing = 1827.1
spacing = 2275.97 # e32p31
# spacing = 2779.25 # g32p32
confstr = "e32p31"
pathstr = f"./data/{confstr}/correlationp2p"

pion_mass_idx = [-1, 0, -6, -7, -13, -12, 17, 16]
rhon_mass_idx = [-3, -2, -5, -9, -15, -8, -11, -14, -10, -19, -18, -4]

p = np.zeros((8, 3))
r = np.zeros((12, 3))

for c in range(20):
    mp, mpe, mm, mme, chidof = fitcorrelation(c, pathstr, nt, t_start, latticespacing=spacing)
    for i in range(len(pion_mass_idx)):
        if c == abs(pion_mass_idx[i]):
            if pion_mass_idx[i] <= 0:
                p[i, 0] = mm
                p[i, 1] = mme
                p[i, 2] = chidof
            else:
                p[i, 0] = mp
                p[i, 1] = mpe
                p[i, 2] = chidof
    for i in range(len(rhon_mass_idx)):
        if c == abs(rhon_mass_idx[i]):
            if rhon_mass_idx[i] <= 0:
                r[i, 0] = mm
                r[i, 1] = mme
                r[i, 2] = chidof
            else:
                r[i, 0] = mp
                r[i, 1] = mpe
                r[i, 2] = chidof

pstr = []
rstr = []
for i in range(8):
    pstr.append(f"${p[i, 0]:.4f}({p[i, 1]*10000:.0f})$ & ${p[i, 0]*spacing:.0f}({p[i, 1]*spacing:.0f})$ & ${p[i, 2]:.2f}$")

for i in range(12):
    rstr.append(f"${r[i, 0]:.4f}({r[i, 1]*10000:.0f})$ & ${r[i, 0]*spacing:.0f}({r[i, 1]*spacing:.0f})$ & ${r[i, 2]:.2f}$")

formatstr = f"""
$1-$ & $\\gamma_5 \\otimes \\tau _5$($\\pi$PS) & {pstr[0]} & $3-$ & $\\gamma_k\\otimes \\tau_k$(VT) & {rstr[0]}\\\\
$0-$ & $\\gamma _4\\gamma _5\\otimes \\tau _4\\tau _5$($\\tilde{{\\pi}}$SC) & {pstr[1]} & $2-$ & $\\gamma_k\\gamma_4\\otimes \\tau_k\\tau_4$(PV) & {rstr[1]}\\\\
$6+$ & $\\gamma_5 \\otimes \\tau _k\\tau _5$($\pi_3$) & {pstr[2]} & $5+$ & $\\gamma_k\\gamma _4\\otimes \\tau_4$($\\rho_3^A$) & {rstr[2]}\\\\
$7+$ & $\\gamma _4\\gamma_5 \\otimes \\tau _l\\tau _m$($\\tilde{{\\pi}}_3$) & {pstr[3]} & $9+$ & $\\gamma_m\\otimes \\tau _l\\tau_m$($\\rho_6^C$) & {rstr[3]}\\\\
$13-$ & $\\gamma_5 \\otimes \\tau _m\\tau _4$ & {pstr[4]} & $15-$ & $\\gamma_l\\gamma _4\\otimes \\tau_k\\tau _4$($\\rho_6^B$) & {rstr[4]}\\\\
$12-$ & $\\gamma _4\\gamma_5 \\otimes \\tau _m$ & {pstr[5]} & $8+$ & $\\gamma_m\\gamma _4\\otimes \\tau _k\\tau_5$($\\rho_6^D$) & {rstr[5]}\\\\
$17+$ & $\\gamma_5 \\otimes \\tau _4$ & {pstr[6]} & $11-$ & $\\gamma_m\\gamma _4\\otimes \\tau_5$ & {rstr[6]}\\\\
$16+$ & $\\gamma _4\\gamma_5 \\otimes \mathbb{{I}}$($\\eta$) & {pstr[7]} & $14-$ & $\\gamma_l\\otimes \\tau_k$($\\rho_6^A$) & {rstr[7]}\\\\
 & & & & & $10-$ & $\\gamma_m\\otimes \\tau _4\\tau_5$ & {rstr[8]}\\\\
 & & & & & $19+$ & $\\gamma_k\\otimes \\tau_l\\tau _m$($\\rho_3^C$) & {rstr[9]}\\\\
 & & & & & $18+$ & $\\gamma_k\\gamma _4\\otimes \\tau_k\\tau _5$($\\rho_3^D$) & {rstr[10]}\\\\
 & & & & & $4+$ & $\\gamma_k\\otimes \mathbb{{I}}$($\omega$) & {rstr[11]}\\\\
"""

print(formatstr)

xlst1 = np.array([0.95 + 0.02 * i for i in range(8)])
xlst2 = np.array([1.93 + 0.02 * i for i in range(12)])

plt.errorbar(xlst1, np.transpose(p)[0]*spacing, yerr=np.transpose(p)[1]*spacing, fmt='o')
plt.errorbar(xlst2, np.transpose(r)[0]*spacing, yerr=np.transpose(r)[1]*spacing, fmt='d')
# plt.locator_params(axis='y', nbins=10)
plt.gca().yaxis.set_major_locator(LinearLocator(numticks=20))
plt.grid(True, axis='y')
plt.ylabel("m (MeV)")
plt.xticks([1, 2], ['$m_{\\pi}$', '$m_{\\rho}$'])
# plt.savefig(f"m{confstr}.pdf")
plt.show()
