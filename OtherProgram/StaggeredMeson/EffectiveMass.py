import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from MesonStructures import all_signs, third_sign
from UsefulFunctions import fitcorrelation, fitcorrelation_with_half

def csymtry(c):
    sizeofc = np.shape(c)
    retv = np.zeros((sizeofc[0], sizeofc[1] // 2 + 1))
    for i in range(1, sizeofc[1] // 2):
        retv[:, i] = 0.5 * (c[:, i] + c[:, sizeofc[1] - i])
    retv[:, 0] = c[:, 0]
    retv[:, sizeofc[1] // 2] = c[:, sizeofc[1] // 2]
    return retv

#########################################################
"""
R(t)=(C(t)+C(t+1))/(C(t-1)+C(t))
with t = 1, ..., (Nt / 2 - 1)

Solve m(t) as

R(t) = (cosh(m(t-Nt/2+1/2)))/(cosh(t-Nt/2-1/2))
"""
def R(c, nt):
    l = 1 + nt // 2
    return (c[1:(l-1)] + c[2:l]) / (c[1:(l-1)] + c[0:(l-2)])

def solve_staggered_mass(R_t, t, Nt):
    """
    Solves for the effective mass 'm' given the ratio R(t).

    Parameters:
    R_t : float : The ratio (C(t) + C(t+1)) / (C(t-1) + C(t))
    t   : int   : The time slice
    Nt  : int   : The total temporal extent of the lattice

    Returns:
    m   : float : The effective mass
    """

    # Define the transcendental equation f(m) = 0
    # We use the midpoint shift for staggered fermions: t - Nt/2
    # R_t = cosh(m * (t - Nt/2 + 0.5)) / cosh(m * (t - Nt/2 - 0.5))

    def f(m):
        if m == 0: return 1.0 - R_t  # Limit as m -> 0
        term1 = np.cosh(m * (t - Nt / 2 + 0.5))
        term2 = np.cosh(m * (t - Nt / 2 - 0.5))
        return term1 - R_t * term2

    # Initial guess using the log approximation to help define brackets
    # m ~ -log(R_t) but we need a safe range for brentq
    try:
        # We search for m in a physically reasonable range [1e-5, 5.0]
        # Most lattice masses fall well within this range.
        m_sol = brentq(f, 1e-8, 5.0)
        return m_sol
    except ValueError:
        # If the ratio is physically impossible or noisy, brentq may fail
        return np.nan

def oneraw_sym_gemini1(c, nt):
    rt = R(c, nt)
    masslst = []
    for i in range(len(rt)):
        masslst.append(solve_staggered_mass(rt[i], i + 1, nt))
    return np.array(masslst)

#########################################################
"""
Cp=C(t)+C(t+1)
Cm=(-1)^t(C(t)-C(t+1))

use one of them

"""

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

def cp_onerow(csym):
    ret = []
    for i in range(len(csym) - 1):
        ret.append(csym[i] + csym[i + 1])
    return np.array(ret)

def cm_onerow(csym):
    ret = []
    for i in range(len(csym) - 1):
        ret.append((csym[i] - csym[i + 1]) * ((-1)**i))
    return np.array(ret)

def cpm_to_mass_1(cpm):
    r = (cpm[2:] + cpm[:-2]) / (2.0 * cpm[1:-1])
    return np.arccosh(r)

def solve_staggered_mass2(R_t, t, Nt):
    # R_t = cosh(m * (t - Nt/2 + 0.5)) / cosh(m * (t - Nt/2 + 1.5))
    def f(m):
        if m == 0: return 1.0 - R_t  # Limit as m -> 0
        term1 = np.cosh(m * (t - Nt / 2 + 0.5))
        term2 = np.cosh(m * (t - Nt / 2 + 1.5))
        return term1 - R_t * term2

    # Initial guess using the log approximation to help define brackets
    # m ~ -log(R_t) but we need a safe range for brentq
    try:
        # We search for m in a physically reasonable range [1e-5, 5.0]
        # Most lattice masses fall well within this range.
        m_sol = brentq(f, 1e-8, 5.0)
        return m_sol
    except ValueError:
        # If the ratio is physically impossible or noisy, brentq may fail
        return np.nan

def cpm_to_mass_2(cpm, nt):
    """
    R(t)=W(t)/W(t+1)
    solve R(t)= (cosh(m(t-Nt/2+1/2)))/(cosh(t-Nt/2+3/2))
    """
    r = cpm[:-1]/cpm[1:]
    masslst = []
    for i in range(len(r)):
        masslst.append(solve_staggered_mass2(r[i], i + 1, nt))
    return np.array(masslst)

def onerow_mass_p(csyn, nt):
    cpm = cp_onerow(csyn)
    r = (cpm[2:] + cpm[:-2]) / (2.0 * cpm[1:-1])
    return np.arccosh(r)

def onerow_mass_m(csyn, nt):
    cpm = cm_onerow(csyn)
    r = (cpm[2:] + cpm[:-2]) / (2.0 * cpm[1:-1])
    return np.arccosh(r)

def onerow_mass_p2(csyn, nt):
    cpm = cp_onerow(csyn)
    return cpm_to_mass_2(cpm, nt)

def onerow_mass_m2(csyn, nt):
    cpm = cm_onerow(csyn)
    return cpm_to_mass_2(cpm, nt)

#########################################################

#########################################################
"""
average test:

jacknife:
"""
def jackknife_meff_func(csym, nt, onerowmassfunc):
    Ncfg = csym.shape[0]

    # full ensemble average
    Cmean = np.mean(csym, axis=0)
    theta_full = onerowmassfunc(Cmean, nt)

    # jackknife samples
    theta_jk = []

    for k in range(Ncfg):
        C_jk = (np.sum(csym, axis=0) - csym[k]) / (Ncfg - 1)
        theta_jk.append(onerowmassfunc(C_jk, nt))
    theta_jk = np.array(theta_jk)
    tvals = np.arange(1, theta_jk.shape[1] + 1)

    # jackknife mean
    theta_bar = np.mean(theta_jk, axis=0)

    # bias correction
    theta_bias_corrected = Ncfg * theta_full - (Ncfg - 1) * theta_bar

    # jackknife variance
    var = (Ncfg - 1) / Ncfg * np.sum((theta_jk - theta_bar)**2, axis=0)

    err = np.sqrt(var)

    return theta_bias_corrected, err, tvals

#########################################################

nt = 64
pathstr = "./data/g32p32/correlationp2p"

all_signs_lst = all_signs()

def jackknife_meff(c):
    """
    Compute effective mass and jackknife error.

    Parameters
    ----------
    c : ndarray
        Correlator array with shape (Ncfg, Nx),
        where Nx corresponds to 0 <= x < Lx//2.

    Returns
    -------
    meff_mean : ndarray
        Bias-corrected jackknife mean of effective mass.
    meff_err : ndarray
        Jackknife error.
    xvals : ndarray
        Positions where m_eff is defined (1 .. Nx-2).
    """

    Ncfg, Nx = c.shape

    # valid x for arccosh estimator
    xvals = np.arange(1, Nx-1)

    # full ensemble average
    Cmean = np.mean(c, axis=0)

    def compute_meff(C):
        r = np.abs((C[2:] + C[:-2]) / (2.0 * C[1:-1]))
        return np.arccosh(r)

    theta_full = compute_meff(Cmean)

    # jackknife samples
    theta_jk = np.zeros((Ncfg, Nx-2))

    for k in range(Ncfg):
        C_jk = (np.sum(c, axis=0) - c[k]) / (Ncfg - 1)
        theta_jk[k] = compute_meff(C_jk)

    # jackknife mean
    theta_bar = np.mean(theta_jk, axis=0)

    # bias correction
    theta_bias_corrected = Ncfg * theta_full - (Ncfg - 1) * theta_bar

    # jackknife variance
    var = (Ncfg - 1) / Ncfg * np.sum((theta_jk - theta_bar)**2, axis=0)

    err = np.sqrt(var)

    return theta_bias_corrected, err, xvals

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

pmtable = [True, False, True, False, False,
           True, False, True, True, False,
           False, True, True, False, False,
           True, True, False, True, False]

drawchannel = [0, 0, 1, 1, 1,
               1, 0, 0, 1, 1,
               1, 1, 0, 0, 1,
               1, 0, 0, 1, 1]

all_third_channel_phase = third_sign()

for channel in range(0, 20):
    correlation_func = np.load(f"{pathstr}_{channel}_0.npy")
    for i in range(1, len(all_signs_lst[channel])):
        correlation_func = correlation_func + np.load(f"{pathstr}_{channel}_{i}.npy")
    corr_sym = csymtry(np.real(correlation_func))
    # testc = np.mean(corr_sym, axis=0)
    if pmtable[channel]:
        y, err, x = jackknife_meff_func(corr_sym, nt, onerow_mass_p2)
    else:
        y, err, x = jackknife_meff_func(corr_sym, nt, onerow_mass_m2)
    if 0 == drawchannel[channel]:
        axes[0].errorbar(x, y, yerr=err, fmt=':o', markersize=2, capsize=2, linewidth=0.2)
    else:
        axes[1].errorbar(x, y, yerr=err, fmt=':o', markersize=2, capsize=2, linewidth=0.2)

    # print(testm)

axes[0].set_xlabel('$x/a$')
axes[0].set_ylabel('$m_{eff}$')
axes[0].set_xlim(0, nt // 2 - 1)
axes[1].set_xlabel('$x/a$')
axes[1].set_ylabel('$m_{eff}$')
axes[1].set_xlim(0, nt // 2 - 1)
plt.show()



