import numpy as np
from matplotlib import pyplot as plt

from CorrelatedFitLMFIT import correlated_fit_with_flexible_bounds_lmfit
from MesonStructures import all_signs


def jackknife_statistic(data, axis=0):
    """
    calculate jackknife statistic
    :param data:
    :param axis:
    :return:
    """
    n = data.shape[axis]
    indices = np.arange(n)
    jackknife_samples = []
    for i in range(n):
        mask = indices != i
        sample = np.take(data, indices[mask], axis=axis)
        jackknife_samples.append(np.mean(sample, axis=axis))
    jackknife_samples = np.array(jackknife_samples)
    original_stat = np.mean(data, axis=axis)
    jackknife_estimate = np.mean(jackknife_samples, axis=0)
    bias = (n - 1) * (jackknife_estimate - original_stat)
    variance = (n - 1) / n * np.sum((jackknife_samples - jackknife_estimate) ** 2, axis=0)
    return {
        'original': original_stat,
        'jackknife_estimate': jackknife_estimate,
        'bias': bias,
        'variance': variance,
        'std_error': np.sqrt(variance),
        'samples': jackknife_samples
    }

def fit_func_draw(nt, t, Ap, mp, Am, mm):
    return Ap * (np.exp(-mp * t) + np.exp(-mp * (nt - t))) + np.cos(np.pi * t) * Am * (np.exp(-mm * t) + np.exp(-mm * (nt - t)))

def draw_fit(ct, ap, mp, am, mm):
    """

    :param ct: ct[conf, nt]
    :param ap: fitted
    :param mp:
    :param am:
    :param mm:
    :return:
    """
    xx_data = []
    yy_data = []
    nt = len(ct[0])
    jk_res = jackknife_statistic(np.real(ct))
    x_data = np.array([t for t in range(nt)])
    y_data = jk_res['jackknife_estimate']
    for t in range(nt):
        xx_data.append(t)
        yy_data.append(fit_func_draw(nt, t, ap, mp, am, mm))
    plt.plot(x_data[1:], y_data[1:])
    plt.plot(xx_data[1:], yy_data[1:], 'o')
    plt.show()

def fitcorrelation(type: int, path: str, nt: int, t_start: int, t_end: int = -1, latticespacing = 1.0, method: str = 'least_squares'):
    """

    """
    all_signs_lst = all_signs()
    all_type = len(all_signs_lst[type])

    def my_func(t, args):
        Ap = args[0]
        mp = args[1]
        Am = args[2]
        mm = args[3]
        return Ap * (np.exp(-mp * t) + np.exp(-mp * (nt - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (nt - t)))

    correlation_func = np.load(f"{path}_{type}_0.npy")
    for i in range(1, all_type):
        correlation_func = correlation_func + np.load(f"{path}_{type}_{i}.npy")

    correlation_func = np.real(correlation_func)
    # pick a separation of t
    t_value = np.array([k for k in range(nt)])
    if t_end > 0:
        t_stride = t_end - t_start
        pickidx = [k + t_start + 1 for k in range(t_stride)] + [nt - t_stride - t_start + k for k in range(t_stride)]
    else:
        pickidx = [k + t_start + 1 for k in range(nt - 1 - 2 * t_start)]
    x_pick = t_value[pickidx]
    y_pick = correlation_func[:, pickidx]
    initial_guess = np.array([np.mean(correlation_func), 0.1, np.std(correlation_func), 0.1])
    upper_bounds = [np.inf, 2.0, np.inf, 2.0]
    lower_bounds = [-np.inf, 0.05, -np.inf, 0.05]

    """
        method : str
            lmfit 优化方法，应该与scipy的curve_fit对应：
            - 'least_squares': 最小二乘法（默认，对应'trf'/'dogbox'）
            - 'leastsq': Levenberg-Marquardt（对应'lm'）
            - 'lbfgsb': L-BFGS-B 算法
            - 'differential_evolution': 差分进化法（全局优化）
            - 'nelder': Nelder-Mead 单纯形法
    """

    fitres = correlated_fit_with_flexible_bounds_lmfit(
            y_pick,
            x_pick,
            my_func,
            initial_guess,
            param_bounds=(lower_bounds, upper_bounds),
            method=method
    )
    print(f"channel:{type}, m+ {fitres['params'][1]} +- {fitres['params_err'][1]} ({fitres['params'][1]*latticespacing} +- {fitres['params_err'][1]*latticespacing}) and m- = {fitres['params'][3]} +- {fitres['params_err'][3]} ({fitres['params'][3]*latticespacing} +- {fitres['params_err'][3]*latticespacing}), chi2/dof = {fitres['chi2_dof']}")
    return fitres['params'][1], fitres['params_err'][1], fitres['params'][3], fitres['params_err'][3], fitres['chi2_dof']

def fitcorrelation_with_half(type: int, path: str, nt: int, t_start: int, method: str = 'least_squares'):
    all_signs_lst = all_signs()
    all_type = len(all_signs_lst[type])

    def my_func(t, args):
        Ap = args[0]
        mp = args[1]
        Am = args[2]
        mm = args[3]
        return Ap * (np.exp(-mp * t) + np.exp(-mp * (nt - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (nt - t)))

    correlation_func = np.load(f"{path}_{type}_0.npy")
    for i in range(1, all_type):
        correlation_func = correlation_func + np.load(f"{path}_{type}_{i}.npy")

    correlation_func = np.real(correlation_func)

    con_count = np.shape(correlation_func)[0]
    nt_count = np.shape(correlation_func)[1]
    data_half = np.zeros((con_count, nt_count // 2))
    for i in range(1, nt // 2 + 1):
        if i == nt // 2:
            data_half[:, i - 1] = correlation_func[:, i]
        else:
            data_half[:, i - 1] = (correlation_func[:, i] + correlation_func[:, nt - i]) / 2
    # data_half_mean = np.mean(data_half, axis=0)
    # print(data_half_mean)
    # plt.plot(data_half_mean)
    # plt.show()
    t_value = np.array([k for k in range(nt)])
    x_pick = t_value[t_start:(t_start+4)]
    y_pick = data_half[:, t_start:(t_start+4)]

    initial_guess = np.array([np.mean(correlation_func), 0.1, np.std(correlation_func), 0.1])
    upper_bounds = [np.inf, 2.0, np.inf, 2.0]
    lower_bounds = [-np.inf, 0.05, -np.inf, 0.05]
    """
        method : str
            lmfit 优化方法，应该与scipy的curve_fit对应：
            - 'least_squares': 最小二乘法（默认，对应'trf'/'dogbox'）
            - 'leastsq': Levenberg-Marquardt（对应'lm'）
            - 'lbfgsb': L-BFGS-B 算法
            - 'differential_evolution': 差分进化法（全局优化）
            - 'nelder': Nelder-Mead 单纯形法
    """
    fitres = correlated_fit_with_flexible_bounds_lmfit(
            y_pick,
            x_pick,
            my_func,
            initial_guess,
            param_bounds=(lower_bounds, upper_bounds),
            method=method
    )
    print(f"channel:{type}, m+ {fitres['params'][1]} +- {fitres['params_err'][1]} and m- = {fitres['params'][3]} +- {fitres['params_err'][3]}, chi2/dof = {fitres['chi2_dof']}")
    return fitres['params'][1], fitres['params_err'][1], fitres['params'][3], fitres['params_err'][3]