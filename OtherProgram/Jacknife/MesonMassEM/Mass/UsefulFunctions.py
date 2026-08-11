import numpy as np
from scipy.optimize import minimize

from MesonMassEM.Mass.CorrelatedFit import correlated_fit_with_flexible_bounds_lmfit


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
        r = (C[2:] + C[:-2]) / (2.0 * C[1:-1])
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

def jackknife_meff2(c):
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
        r = np.abs(C[:-2] / C[2:])
        return 0.5 * np.log(r)

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

def staggered_correlated_fit(C, Lx, x_min, fit_osc=0):
    """

    :param C:
    :param Lx:
    :param x_min:
    :param fit_osc: 0, both, 1, only positive, 2, only negative
    :return:
    """

    Ncfg, Nx = C.shape
    x_max = Lx // 2
    x = np.arange(x_min, x_max + 1)
    idx = slice(x_min, x_max + 1)

    # ensemble mean
    Cmean = np.mean(C[:, idx], axis=0)

    # covariance
    delta = C[:, idx] - Cmean
    cov = delta.T @ delta / (Ncfg - 1)

    # SVD inverse
    U, s, Vt = np.linalg.svd(cov)
    cutoff = s.max() * 1e-8
    s_inv = np.array([1/si if si > cutoff else 0 for si in s])
    cov_inv = (Vt.T * s_inv) @ U.T

    xc = x - Lx/2

    # model
    def model(p):

        if 0 == fit_osc:

            Aplus, muplus, Aminus, muminus = p
            mplus = np.exp(muplus)
            mminus = np.exp(muminus)

            return (
                Aplus * np.cosh(mplus * xc)
                + ((-1)**x) * Aminus * np.cosh(mminus * xc)
            )

        elif 1 == fit_osc:

            Aplus, muplus = p
            mplus = np.exp(muplus)

            return Aplus * np.cosh(mplus * xc)
        else:
            Aminus, muminus = p
            mminus = np.exp(muminus)

            return ((-1)**x) * Aminus * np.cosh(mminus * xc)

    # chi2
    def chi2(p, Cdata):

        y = model(p)

        if not np.all(np.isfinite(y)):
            return 1e100

        r = Cdata - y
        return r @ cov_inv @ r

    # effective mass guess
    try:
        meff = abs(np.arccosh((Cmean[-2] + Cmean[-4])/(2*Cmean[-3])))
    except:
        meff = 0.6

    if not np.isfinite(meff) or meff <= 0:
        meff = 0.6

    # initial parameters
    if 0 == fit_osc:

        p0 = np.array([
            Cmean[-1],
            np.log(meff),
            0.1*Cmean[-1],
            np.log(1.5*meff)
        ])

        bounds = [
            (None,None),
            (-5,np.log(5)),
            (None,None),
            (-5,np.log(5))
        ]

    else:

        p0 = np.array([
            Cmean[-1],
            np.log(meff)
        ])

        bounds = [
            (None,None),
            (-5,np.log(5))
        ]

    # central fit
    res = minimize(
        chi2,
        p0,
        args=(Cmean,),
        bounds=bounds,
        method="L-BFGS-B"
    )

    if not res.success:
        raise RuntimeError("Fit failed")

    p_fit = res.x

    if 0 == fit_osc:

        Aplus = p_fit[0]
        mplus = np.exp(p_fit[1])
        Aminus = p_fit[2]
        mminus = np.exp(p_fit[3])

        dof = len(x) - 4

    elif 1 == fit_osc:

        Aplus = p_fit[0]
        mplus = np.exp(p_fit[1])
        Aminus = np.nan
        mminus = np.nan

        dof = len(x) - 2
    else:
        Aplus = np.nan
        mplus = np.nan
        Aminus = p_fit[0]
        mminus = np.exp(p_fit[1])

        dof = len(x) - 2

    chi2_dof = chi2(p_fit, Cmean) / dof

    # jackknife
    sumC = np.sum(C[:, idx], axis=0)
    params = []

    for k in range(Ncfg):

        Cjk = (sumC - C[k, idx]) / (Ncfg - 1)

        res = minimize(
            chi2,
            p_fit,
            args=(Cjk,),
            bounds=bounds,
            method="L-BFGS-B"
        )

        if res.success:
            params.append(res.x)

    params = np.array(params)

    if 0 == fit_osc:

        mplus_jk = np.exp(params[:,1])
        Aplus_jk = params[:,0]
        mminus_jk = np.exp(params[:,3])
        Aminus_jk = params[:,2]

    elif 1 == fit_osc:

        mplus_jk = np.exp(params[:,1])
        Aplus_jk = params[:,0]
    else:
        mminus_jk = np.exp(params[:,1])
        Aminus_jk = params[:,0]



    if 0 == fit_osc:
        mplus_mean = np.mean(mplus_jk)
        Aplus_mean = np.mean(Aplus_jk)

        dmplus = np.sqrt((Ncfg - 1) * np.var(mplus_jk))
        dAplus = np.sqrt((Ncfg - 1) * np.var(Aplus_jk))

        mminus_mean = np.mean(mminus_jk)
        Aminus_mean = np.mean(Aminus_jk)

        dmminus = np.sqrt((Ncfg-1)*np.var(mminus_jk))
        dAminus = np.sqrt((Ncfg-1)*np.var(Aminus_jk))

    elif 1 == fit_osc:
        mplus_mean = np.mean(mplus_jk)
        Aplus_mean = np.mean(Aplus_jk)

        dmplus = np.sqrt((Ncfg - 1) * np.var(mplus_jk))
        dAplus = np.sqrt((Ncfg - 1) * np.var(Aplus_jk))

        mminus_mean = np.nan
        Aminus_mean = np.nan
        dmminus = np.nan
        dAminus = np.nan
    else:

        mplus_mean = np.nan
        Aplus_mean = np.nan
        dmplus = np.nan
        dAplus = np.nan

        mminus_mean = np.mean(mminus_jk)
        Aminus_mean = np.mean(Aminus_jk)

        dmminus = np.sqrt((Ncfg-1)*np.var(mminus_jk))
        dAminus = np.sqrt((Ncfg-1)*np.var(Aminus_jk))


    return (
        mplus_mean, Aplus_mean,
        mminus_mean, Aminus_mean,
        dmplus, dAplus,
        dmminus, dAminus,
        chi2_dof
    )