import numpy as np
from scipy.optimize import minimize

def positive_channel(t, params):
    T = 24
    return params[1] * (np.exp(-params[0] * t) + np.exp(-params[0] * (T - t))) * 0.5

def positive_channel_E(t, params, k):
    T = 24
    m = params[0]
    return params[1] * (np.exp(-m * t) + np.exp(-m * (T - t))) * 0.5 * (1 + params[2] * np.cos(2 * k * np.pi * t / T))

def positive_channel_E2(t, params, k):
    T = 24
    m = params[0] * np.sqrt(1 + params[2] * np.cos(2 * k * np.pi * t / T) + params[3] * np.cos(4 * k * np.pi * t / T))
    return params[1] * (np.exp(-m * t) + np.exp(-m * (T - t))) * 0.5

def positive_channel_1(t, params):
    return positive_channel_E(t, params, 1)

def positive_channel_2(t, params):
    return positive_channel_E(t, params, 2)

def positive_channel_3(t, params):
    return positive_channel_E(t, params, 3)

def positive_channel_4(t, params):
    return positive_channel_E(t, params, 4)

def positive_channel_5(t, params):
    return positive_channel_E(t, params, 5)

def positive_channel_6(t, params):
    return positive_channel_E(t, params, 6)

def positive_channel_7(t, params):
    return positive_channel_E(t, params, 7)

def positive_channel_8(t, params):
    return positive_channel_E(t, params, 8)

def positive_channel_10(t, params):
    return positive_channel_E(t, params, 10)

def positive_channel_12(t, params):
    return positive_channel_E(t, params, 12)

def positive_channel_14(t, params):
    return positive_channel_E(t, params, 14)

def positive_channelf_1(t, params):
    return positive_channel_E2(t, params, 1)

def positive_channelf_2(t, params):
    return positive_channel_E2(t, params, 2)

def positive_channelf_3(t, params):
    return positive_channel_E2(t, params, 3)

def positive_channelf_4(t, params):
    return positive_channel_E2(t, params, 4)

def positive_channelf_5(t, params):
    return positive_channel_E2(t, params, 5)

def positive_channelf_6(t, params):
    return positive_channel_E2(t, params, 6)

def positive_channelf_7(t, params):
    return positive_channel_E2(t, params, 7)

def negative_channel(t, params):
    T = 24
    return ((-1)**t) * params[1] * (np.exp(-params[0] * t) + np.exp(-params[0] * (T - t))) * 0.5

def negative_channel_draw(t, params):
    T = 24
    return np.cos(np.pi * t) * params[1] * (np.exp(-params[0] * t) + np.exp(-params[0] * (T - t))) * 0.5

def chiral_fit_1(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(2 * np.pi * t / 24)

def chiral_fit_2(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(4 * np.pi * t / 24)

def chiral_fit_3(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(6 * np.pi * t / 24)

def chiral_fit_4(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(8 * np.pi * t / 24)

def chiral_fit_5(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(10 * np.pi * t / 24)

def chiral_fit_6(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(12 * np.pi * t / 24)

def chiral_fit_7(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(14 * np.pi * t / 24)

def chiral_fit_8(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(16 * np.pi * t / 24)

def chiral_fit_10(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(20 * np.pi * t / 24)

def chiral_fit_12(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(24 * np.pi * t / 24)

def chiral_fit_14(t, params):
    a = params[0]
    b = params[1]
    return a + b * np.cos(28 * np.pi * t / 24)

def charge_fit_1(t, params):
    b = params[0]
    return b * np.sin(2 * np.pi * t / 24)

def charge_fit_2(t, params):
    b = params[0]
    return b * np.sin(4 * np.pi * t / 24)

def charge_fit_3(t, params):
    b = params[0]
    return b * np.sin(6 * np.pi * t / 24)

def charge_fit_4(t, params):
    b = params[0]
    return b * np.sin(8 * np.pi * t / 24)

def charge_fit_5(t, params):
    b = params[0]
    return b * np.sin(10 * np.pi * t / 24)

def charge_fit_6(t, params):
    b = params[0]
    return b * np.sin(12 * np.pi * t / 24)

def charge_fit_7(t, params):
    b = params[0]
    return b * np.sin(14 * np.pi * t / 24)

def charge_fit_8(t, params):
    b = params[0]
    return b * np.sin(16 * np.pi * t / 24)

def charge_fit_10(t, params):
    b = params[0]
    return b * np.sin(20 * np.pi * t / 24)

def charge_fit_12(t, params):
    b = params[0]
    return b * np.sin(24 * np.pi * t / 24)

def charge_fit_14(t, params):
    b = params[0]
    return b * np.sin(28 * np.pi * t / 24)

def polyakov_fit(t, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a + b * np.exp(-2j * np.pi * t / 24) + c * np.exp(4j * np.pi * t / 24)

def polyakov_fit2(t, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a + b * np.exp(-4j * np.pi * t / 24) + c * np.exp(8j * np.pi * t / 24)

def polyakov_fit580b(t, params):
    a = params[0]
    b = params[1]
    c = params[2]
    d = params[3]
    e = params[4]
    return (a + d * np.cos(2 * np.pi * t / 24) + e * np.cos(4 * np.pi * t / 24) + b * np.exp(-2j * np.pi * t / 24) + c * np.exp(4j * np.pi * t / 24))

def polyakov_fit580b2(t, params):
    a = params[0]
    b = params[1]
    c = params[2]
    d = params[3]
    e = params[4]
    return (a + d * np.cos(4 * np.pi * t / 24) + e * np.cos(8 * np.pi * t / 24) + b * np.exp(-4j * np.pi * t / 24) + c * np.exp(8j * np.pi * t / 24))

def mp_channel(t, params):
    m0, A0, m1, A1 = params
    T = 24
    return 0.5 * (
        A0 * (np.exp(-m0*t) + np.exp(-m0*(T-t))) +
        ((-1)**t) * A1 * (np.exp(-m1*t) + np.exp(-m1*(T-t)))
    )

def pp_channel(t, params):
    m0, A0, m1, A1 = params
    T = 24
    return 0.5 * (
        A0 * (np.exp(-m0*t) + np.exp(-m0*(T-t))) +
        A1 * (np.exp(-m1*t) + np.exp(-m1*(T-t)))
    )

def correlated_fit(
    c,
    t_vals,
    f,
    p0,
    lower_bounds,
    upper_bounds,
    t_min=None,
    t_max=None,
    svd_cut=None,
    eps=1e-12,
    jackknife=True,
    complex_mode="auto"  # "auto", "real", "imag", "full"
):
    """
    General correlated fit supporting real/complex data and model.

    Parameters
    ----------
    c : array (Ncfg, Nt), real or complex
    t_vals : array (Nt,), real
    f : callable, f(t, params) -> real or complex
    p0 : initial parameters (real)
    lower_bounds, upper_bounds : bounds (np.nan = no bound)
    t_min, t_max : fit window
    svd_cut : eigenvalue cutoff
    eps : regularization
    jackknife : bool
    complex_mode :
        "auto"  -> detect from data
        "real"  -> use Re only
        "imag"  -> use Im only
        "full"  -> use full complex

    Returns
    -------
    params : best-fit parameters
    param_err : jackknife errors (or None)
    chi2_dof : chi^2 / dof
    """

    c = np.asarray(c)
    t_vals = np.asarray(t_vals)

    Ncfg, Nt = c.shape
    k = len(p0)

    # ----------------------------
    # (A) Fit window
    # ----------------------------
    if t_min is not None or t_max is not None:
        mask = np.ones_like(t_vals, dtype=bool)
        if t_min is not None:
            mask &= (t_vals >= t_min)
        if t_max is not None:
            mask &= (t_vals <= t_max)

        t_vals = t_vals[mask]
        c = c[:, mask]
        Nt = len(t_vals)

    # ----------------------------
    # determine complex handling
    # ----------------------------
    if complex_mode == "auto":
        if np.iscomplexobj(c):
            imag_ratio = np.mean(np.abs(c.imag)) / (np.mean(np.abs(c.real)) + 1e-16)
            complex_mode = "real" if imag_ratio < 1e-8 else "full"
        else:
            complex_mode = "real"

    # ----------------------------
    # helper: embed data
    # ----------------------------
    def embed(vec):
        if complex_mode == "real":
            return vec.real
        elif complex_mode == "imag":
            return vec.imag
        elif complex_mode == "full":
            return np.concatenate([vec.real, vec.imag])
        else:
            raise ValueError("Invalid complex_mode")

    # ----------------------------
    # covariance inverse
    # ----------------------------
    def compute_cov_inv(data):
        mean = np.mean(data, axis=0)
        fluctuations = data - mean

        if complex_mode == "full":
            cov = fluctuations.conj().T @ fluctuations / (len(data) - 1)
            cov /= len(data)

            Nt_loc = len(mean)

            cov_real = np.block([
                [cov.real, -cov.imag],
                [cov.imag, cov.real]
            ])
        else:
            data_proj = embed(data)
            mean = np.mean(data_proj, axis=0)
            fluctuations = data_proj - mean

            cov_real = fluctuations.T @ fluctuations / (len(data_proj) - 1)
            cov_real /= len(data_proj)

        cov_real += eps * np.eye(len(cov_real))

        # ---- SVD part ----
        if svd_cut is not None:
            U, s, Vt = np.linalg.svd(cov_real)
            mask = s > svd_cut
            N_eff = np.sum(mask)

            if N_eff == 0:
                raise RuntimeError("All singular values removed by svd_cut")

            cov_inv = (Vt.T[:, mask] * (1.0 / s[mask])) @ U.T[mask, :]
        else:
            N_eff = cov_real.shape[0]
            try:
                cov_inv = np.linalg.inv(cov_real)
            except np.linalg.LinAlgError:
                # cov_inv = np.linalg.pinv(cov_real)

                # ---- 显式 SVD，替代 pinv ----
                U, s, Vt = np.linalg.svd(cov_real)
                # 数值秩判据（等价于 pinv 内部逻辑）
                tol = np.max(s) * max(cov_real.shape) * np.finfo(float).eps
                mask = s > tol
                N_eff = np.sum(mask)
                if N_eff == 0:
                    raise RuntimeError("Covariance matrix is numerically zero")
                cov_inv = (Vt.T[:, mask] * (1.0 / s[mask])) @ U.T[mask, :]

        return mean, cov_inv, N_eff

    # ----------------------------
    # bounds
    # ----------------------------
    bounds = []
    for lo, hi in zip(lower_bounds, upper_bounds):
        lo = None if np.isnan(lo) else lo
        hi = None if np.isnan(hi) else hi
        bounds.append((lo, hi))

    # ----------------------------
    # single fit
    # ----------------------------
    def do_fit(data):
        mean, cov_inv, N_eff = compute_cov_inv(data)

        def chi2(params):
            model = np.array([f(t, params) for t in t_vals])

            # if not np.all(np.isfinite(model)):
            #     return 1e100  # reject bad region

            if complex_mode != "full":
                model = embed(model)
                diff = mean - model
            else:
                diff = mean - model
                diff = np.concatenate([diff.real, diff.imag])

            val = diff @ cov_inv @ diff

            # if not np.isfinite(val):
            #     return 1e100

            return val

        result = minimize(chi2, p0, bounds=bounds)

        if not result.success:
            raise RuntimeError("Fit failed: " + result.message)

        return result.x, result.fun, N_eff

    # ----------------------------
    # central fit
    # ----------------------------
    params, chi2_val, N_eff = do_fit(c)

    dof = N_eff - k
    chi2_dof = chi2_val / dof

    # ----------------------------
    # (B) Jackknife
    # ----------------------------
    param_err = None

    if jackknife:
        jk_params = []

        for i in range(Ncfg):
            mask = np.ones(Ncfg, dtype=bool)
            mask[i] = False
            data_jk = c[mask]

            try:
                p_jk, _, _ = do_fit(data_jk)
                # ---- sanity check ----
                # if not np.all(np.isfinite(p_jk)):
                #     continue

                # optional: reject unreasonable mass values
                # if np.any(np.abs(p_jk) > 10):  # adjust threshold
                #     continue

                jk_params.append(p_jk)
            except RuntimeError:
                continue

        jk_params = np.array(jk_params)

        if len(jk_params) < Ncfg // 2:
            raise RuntimeError("Too many jackknife failures")

        mean_jk = np.mean(jk_params, axis=0)

        param_err = np.sqrt(
            (len(jk_params) - 1) *
            np.mean((jk_params - mean_jk) ** 2, axis=0)
        )

    return params, param_err, chi2_dof

def MPFit_P(c, x_min, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel, [0.5, 1], [0.01, np.nan], [2.0, np.nan], x_min, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], chi2

def MPFit_M(c, x_min, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), negative_channel, [0.5, 1], [0.01, np.nan], [2.0, np.nan], x_min, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], chi2

def MPFit_ED(c, x_min, k, svd_cut=None):
    shapec = np.shape(c)
    if 0 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel, [0.5, 1], [0.01, np.nan],[1.0, np.nan], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], 0, 0, chi2
    elif 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_1, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_2, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_3, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_4, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_5, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_6, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_7, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    return 0,0,0,0,0,0,0

def MPFit_EU(c, x_min, k, svd_cut=None):
    shapec = np.shape(c)
    if 0 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel, [0.5, 1], [0.01, np.nan],[1.0, np.nan], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], 0, 0, chi2
    elif 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_2, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_4, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_6, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_8, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_10, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_12, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel_14, [0.5, 1, 0], [0.01, np.nan, -0.5], [1.0, np.nan, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2
    return 0,0,0,0,0,0,0

def MPFit_EM(c, x_min, k, svd_cut=None):
    shapec = np.shape(c)
    if k < 8:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channel, [0.5, 1], [0.01, np.nan],[1.0, np.nan], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], 0, 0, 0, 0, chi2
    """
    elif 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_1, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_2, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_3, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_4, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_5, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_6, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), positive_channelf_7, [0.5, 1, 0, 0], [0.01, np.nan, -0.5, -0.5], [1.0, np.nan, 0.5, 0.5], x_min, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], chi2
    """
    return 0,0,0,0,0,0,0,0,0

def M_PMFit(c, x_min, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), mp_channel, [0.5, 1, 0.5, 1], [0.01, np.nan, 0.01, np.nan], [5.0, np.nan, 5.0, np.nan], x_min, svd_cut=svd_cut)
    return params[0], param_err[0], params[2], param_err[2], chi2

def M_PPFit(c, x_min, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), pp_channel, [0.5, 1, 0.5, 1], [0.01, np.nan, 0.01, np.nan], [2.0, np.nan, 2.0, np.nan], x_min, svd_cut=svd_cut)
    return params[0], param_err[0], params[2], param_err[2], chi2

def ChiralFit_D(c, k, svd_cut=None):
    shapec = np.shape(c)
    if 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_1, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_2, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_3, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_4, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_5, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_6, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_7, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    return 0,0,0,0,0

def ChiralFit_U(c, k, svd_cut=None):
    shapec = np.shape(c)
    if 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_2, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_4, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_6, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_8, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_10, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_12, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), chiral_fit_14, [-0.1, 0.1], [-0.5, -0.5], [0.5, 0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], params[1], param_err[1], chi2
    return 0,0,0,0,0

def ChargeFit_D(c, k, svd_cut=None):
    shapec = np.shape(c)
    if 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_1, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_2, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_3, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_4, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_5, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_6, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_7, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    return 0,0,0

def ChargeFit_U(c, k, svd_cut=None):
    shapec = np.shape(c)
    if 1 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_2, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 2 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_4, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 3 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_6, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 4 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_8, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 5 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_10, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 6 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_12, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    elif 7 == k:
        params, param_err, chi2 = correlated_fit(c, range(shapec[1]), charge_fit_14, [0.1], [-0.5], [0.5], 0, svd_cut=svd_cut)
        return params[0], param_err[0], chi2
    return 0,0,0

def PolyaFit(c, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), polyakov_fit, [0.1, 0, 0], [-0.1, -0.1, -0.1], [0.3, 0.1, 0.1], 0, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2

def PolyaFit2(c, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), polyakov_fit2, [0.1, 0, 0], [-0.1, -0.1, -0.1], [0.3, 0.1, 0.1], 0, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], chi2

def PolyaFit580(c, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), polyakov_fit580b, [0.1, 0, 0, 0, 0], [0.1, -0.1, -0.1, -0.8, -0.8], [0.3, 0.1, 0.1, 0.8, 0.8], 0, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], params[4], param_err[4], chi2

def PolyaFit5802(c, svd_cut=None):
    shapec = np.shape(c)
    params, param_err, chi2 = correlated_fit(c, range(shapec[1]), polyakov_fit580b2, [0.1, 0, 0, 0, 0], [0.1, -0.1, -0.1, -0.8, -0.8], [0.3, 0.1, 0.1, 0.8, 0.8], 0, svd_cut=svd_cut)
    return params[0], param_err[0], params[1], param_err[1], params[2], param_err[2], params[3], param_err[3], params[4], param_err[4], chi2