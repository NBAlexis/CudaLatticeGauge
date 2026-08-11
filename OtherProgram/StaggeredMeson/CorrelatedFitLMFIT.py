import numpy as np
import lmfit
from scipy.linalg import inv, pinv, cholesky, solve_triangular
import warnings
from typing import Callable, Optional, Union, List, Tuple


def correlated_fit_with_flexible_bounds_lmfit(
        data_matrix: np.ndarray,
        x_values: np.ndarray,
        model_func: Callable,
        initial_params: np.ndarray,
        param_bounds: Optional[
            Union[List[Tuple[Optional[float], Optional[float]]], Tuple[List[float], List[float]]]] = None,
        method: str = 'least_squares',
        compute_correlation: bool = True,
        verbose: bool = False,
        uncertainty_method: str = 'covariance'
) -> dict:
    """
    使用 lmfit 执行关联拟合，支持灵活的边界约束

    完全使用 lmfit 实现，不依赖 scipy.curve_fit

    参数：
    ----------
    data_matrix : numpy.ndarray
        形状为 (n_configs, m_points) 的矩阵
        每行是一个组态，每列是一个测量点y_i
    x_values : numpy.ndarray
        形状为 (m_points,) 的数组，x_i的值
    model_func : callable
        拟合函数 f(x, params)，返回预测值
        params是待拟合参数列表
    initial_params : list or numpy.ndarray
        参数的初始猜测值
    param_bounds : 多种格式可选
        1. None: 无边界
        2. 传统格式: ([min1, min2, ...], [max1, max2, ...]) - 元组包含两个列表
        3. 灵活格式: [(min1, max1), (min2, max2), ...] - 列表包含元组
        其中min和max可以是：
        - float: 具体的边界值
        - None: 无限制（负无穷或正无穷）
        - np.inf或-np.inf: 无穷大或无穷小
    method : str
        lmfit 优化方法：
        - 'least_squares': 最小二乘法（默认，对应scipy的'trf'/'dogbox'）
        - 'leastsq': Levenberg-Marquardt（对应scipy的'lm'）
        - 'lbfgsb': L-BFGS-B 算法
        - 'differential_evolution': 差分进化法（全局优化）
        - 'nelder': Nelder-Mead 单纯形法
    compute_correlation : bool
        是否计算参数间的相关系数矩阵
    verbose : bool
        是否打印详细过程信息
    uncertainty_method : str
        不确定性估计方法：
        - 'covariance': 协方差矩阵方法（默认）
        - 'bootstrap': Bootstrap方法
        - 'mcmc': MCMC方法（需要emcee包）

    返回：
    ----------
    result : dict
        包含以下键：
        - 'params': 拟合参数值
        - 'params_err': 参数误差
        - 'params_bounds': 参数是否在边界上（列表，True表示在边界上）
        - 'bounds_type': 每个参数的边界类型
        - 'chi2': χ²值
        - 'chi2_dof': χ²/自由度
        - 'cov_matrix': 参数的协方差矩阵
        - 'corr_matrix': 参数间的相关系数矩阵（如果compute_correlation=True）
        - 'data_cov': 数据的协方差矩阵
        - 'success': 是否成功
        - 'message': 优化器消息
    """

    # 输入检查
    n_configs, m_points = data_matrix.shape
    if len(x_values) != m_points:
        raise ValueError(f"x_values长度({len(x_values)})与数据矩阵列数({m_points})不匹配")

    initial_params = np.array(initial_params, dtype=float)
    n_params = len(initial_params)

    # 处理边界 - 支持多种格式（与原始函数完全相同）
    bounds_info = None
    bounds_type = None
    has_bounds = False

    # 首先将 param_bounds 转换为统一格式：列表的列表，每个元素是 [min, max]
    processed_bounds = []

    if param_bounds is not None:
        # 检查边界格式
        if isinstance(param_bounds, tuple) and len(param_bounds) == 2:
            # 传统格式: ([min1, min2, ...], [max1, max2, ...])
            lower_list, upper_list = param_bounds

            if len(lower_list) != n_params or len(upper_list) != n_params:
                raise ValueError(
                    f"传统边界格式: 下界列表长度({len(lower_list)})和上界列表长度({len(upper_list)})必须与参数个数({n_params})相同")

            # 转换为灵活格式
            for i in range(n_params):
                min_val = lower_list[i]
                max_val = upper_list[i]

                # 将np.inf/np.nan转换为None
                if min_val is None or np.isnan(min_val):
                    min_val = None
                elif np.isneginf(min_val):
                    min_val = None
                elif np.isinf(min_val) and min_val < 0:
                    min_val = None

                if max_val is None or np.isnan(max_val):
                    max_val = None
                elif np.isposinf(max_val):
                    max_val = None
                elif np.isinf(max_val) and max_val > 0:
                    max_val = None

                processed_bounds.append((min_val, max_val))

        elif isinstance(param_bounds, list):
            # 灵活格式: [(min1, max1), (min2, max2), ...]
            if len(param_bounds) != n_params:
                raise ValueError(f"灵活边界格式: 边界列表长度({len(param_bounds)})必须与参数个数({n_params})相同")

            # 确保每个元素都是元组
            for i, bound in enumerate(param_bounds):
                if isinstance(bound, tuple) and len(bound) == 2:
                    min_val, max_val = bound

                    # 处理边界值
                    if min_val is None or np.isnan(min_val):
                        min_val = None
                    elif np.isneginf(min_val):
                        min_val = None
                    elif np.isinf(min_val) and min_val < 0:
                        min_val = None

                    if max_val is None or np.isnan(max_val):
                        max_val = None
                    elif np.isposinf(max_val):
                        max_val = None
                    elif np.isinf(max_val) and max_val > 0:
                        max_val = None

                    processed_bounds.append((min_val, max_val))
                else:
                    raise ValueError(f"边界元素{i}不是有效的元组格式: {bound}")
        else:
            raise ValueError(
                f"不支持的边界格式: {type(param_bounds)}。应为None、元组([min...],[max...])或列表[(min,max),...]")

        # 记录边界信息
        bounds_info = []
        bounds_type = []

        for i, (min_val, max_val) in enumerate(processed_bounds):
            # 处理min_val
            if min_val is None:
                lb = -np.inf
            elif np.isinf(min_val) and min_val < 0:
                lb = -np.inf
            else:
                lb = float(min_val)
                if not np.isinf(lb):
                    has_bounds = True

            # 处理max_val
            if max_val is None:
                ub = np.inf
            elif np.isinf(max_val) and max_val > 0:
                ub = np.inf
            else:
                ub = float(max_val)
                if not np.isinf(ub):
                    has_bounds = True

            # 记录边界信息
            bounds_info.append((lb, ub))

            # 记录边界类型
            lb_finite = not np.isinf(lb) and lb != -np.inf
            ub_finite = not np.isinf(ub) and ub != np.inf

            if not lb_finite and not ub_finite:
                bounds_type.append("无限制")
            elif lb_finite and ub_finite:
                bounds_type.append("双边限制")
            elif lb_finite:
                bounds_type.append("下界限制")
            else:
                bounds_type.append("上界限制")

        # 检查初始值是否在边界内
        for i, (init, (lb, ub)) in enumerate(zip(initial_params, bounds_info)):
            if init < lb or init > ub:
                warnings.warn(f"参数{i}的初始值{init}不在边界[{lb}, {ub}]内")
                # 调整到边界内（如果需要）
                if init < lb:
                    if lb == -np.inf:
                        initial_params[i] = init
                    else:
                        initial_params[i] = lb + 0.01 * (ub - lb) if ub != np.inf else lb + 0.01
                else:  # init > ub
                    if ub == np.inf:
                        initial_params[i] = init
                    else:
                        initial_params[i] = ub - 0.01 * (ub - lb) if lb != -np.inf else ub - 0.01
    else:
        # 如果没有边界，创建默认边界
        processed_bounds = [(None, None) for _ in range(n_params)]
        bounds_info = [(-np.inf, np.inf) for _ in range(n_params)]
        bounds_type = ["无限制" for _ in range(n_params)]

    if verbose:
        print(f"数据统计: {n_configs}个组态, {m_points}个测量点, {n_params}个参数")
        print("参数边界类型:")
        for i, (param, btype) in enumerate(zip(initial_params, bounds_type)):
            min_val, max_val = bounds_info[i] if bounds_info else (None, None)
            min_str = "-∞" if min_val == -np.inf else f"{min_val:.4f}"
            max_str = "∞" if max_val == np.inf else f"{max_val:.4f}"
            print(f"  参数{i}: {btype}, 初始值={param:.4f}, 边界=[{min_str}, {max_str}]")

    # 1. 计算平均值和协方差矩阵
    Y_mean = np.mean(data_matrix, axis=0)

    # 关键修正：正确计算平均值的协方差矩阵
    # 平均值的协方差 = 样本协方差 / 组态数
    if n_configs > 1:
        # 计算原始数据的样本协方差（无偏估计）
        sample_cov = np.cov(data_matrix, rowvar=False, ddof=1)
        # 平均值的协方差
        mean_cov = sample_cov / n_configs
    else:
        raise ValueError("至少需要2个组态来计算协方差")

    # 诊断信息：计算简单误差（忽略相关性）
    Y_std = np.std(data_matrix, axis=0, ddof=1)
    Y_err_simple = Y_std / np.sqrt(n_configs)  # 平均值的简单标准误差

    if verbose:
        print(f"平均值的简单标准误差范围: [{Y_err_simple.min():.6e}, {Y_err_simple.max():.6e}]")
        print(f"平均值的协方差矩阵对角线均值: {np.mean(np.diag(mean_cov)):.6e}")
        print(f"协方差矩阵条件数: {np.linalg.cond(mean_cov):.6e}")

    # 改进的协方差矩阵求逆，增加数值稳定性
    def safe_inverse(cov_matrix, regularization=1e-12):
        """安全的矩阵求逆，增加数值稳定性"""
        try:
            # 尝试Cholesky分解验证正定性
            L = cholesky(cov_matrix, lower=True, check_finite=False)
            # 使用Cholesky分解求解逆矩阵
            cov_inv = cholesky_solve(L, np.eye(cov_matrix.shape[0]))
            method_used = "cholesky"
            return cov_inv, method_used
        except:
            try:
                # 尝试标准逆
                cov_inv = inv(cov_matrix)
                method_used = "inverse"
                return cov_inv, method_used
            except:
                try:
                    # 添加正则化项
                    n = cov_matrix.shape[0]
                    regularized_cov = cov_matrix + regularization * np.eye(n) * np.trace(cov_matrix) / n
                    cov_inv = inv(regularized_cov)
                    method_used = "regularized_inverse"
                    return cov_inv, method_used
                except:
                    # 最后使用伪逆
                    cov_inv = pinv(cov_matrix, rcond=1e-12)
                    method_used = "pseudo_inverse"
                    return cov_inv, method_used

    def cholesky_solve(L, B):
        """使用Cholesky分解求解线性系统"""
        # 解 L y = B
        y = solve_triangular(L, B, lower=True)
        # 解 L^T x = y
        x = solve_triangular(L.T, y, lower=False)
        return x

    # 计算协方差矩阵的逆
    cov_inv, inv_method = safe_inverse(mean_cov)

    if verbose:
        print(f"协方差矩阵求逆方法: {inv_method}")

    # 2. 创建 lmfit 参数对象
    param_names = [f'p{i}' for i in range(n_params)]
    params = lmfit.Parameters()

    for i in range(n_params):
        name = param_names[i]
        init_val = initial_params[i]
        min_val, max_val = processed_bounds[i]

        # 设置参数，包括边界
        params.add(name, value=init_val, min=min_val, max=max_val)

    # 3. 根据方法类型定义不同的目标函数
    # 将scipy方法名映射到lmfit方法名
    method_mapping = {
        'trf': 'least_squares',
        'dogbox': 'least_squares',
        'lm': 'leastsq',
        'least_squares': 'least_squares',
        'leastsq': 'leastsq',
        'lbfgsb': 'lbfgsb',
        'differential_evolution': 'differential_evolution',
        'nelder': 'nelder'
    }

    lmfit_method = method_mapping.get(method, 'least_squares')

    # 判断方法类型：是需要残差向量的方法，还是需要标量值的方法
    # 需要残差向量的方法：least_squares, leastsq
    # 需要标量值的方法：lbfgsb, differential_evolution, nelder等
    if lmfit_method in ['least_squares', 'leastsq']:
        # 定义返回残差向量的目标函数
        def objective(params):
            """目标函数，返回加权残差向量（用于最小二乘方法）"""
            # 从 params 字典中提取参数值
            param_values = np.array([params[name].value for name in param_names])

            # 计算模型值
            model_vals = model_func(x_values, param_values)

            # 计算残差
            residuals = Y_mean - model_vals

            # 使用协方差矩阵的Cholesky分解的逆来加权残差
            try:
                L = cholesky(mean_cov, lower=True, check_finite=False)
                weighted_residuals = solve_triangular(L, residuals, lower=True)
            except:
                # 如果Cholesky分解失败，使用对角线权重
                weighted_residuals = residuals / Y_err_simple

            return weighted_residuals
    else:
        # 定义返回标量值的目标函数（用于lbfgsb, differential_evolution等方法）
        def objective(params):
            """目标函数，返回标量chi2值（用于标量优化方法）"""
            # 从 params 字典中提取参数值
            param_values = np.array([params[name].value for name in param_names])

            # 计算模型值
            model_vals = model_func(x_values, param_values)

            # 计算残差
            residuals = Y_mean - model_vals

            # 计算 χ² = r^T C^{-1} r
            chi2 = residuals.T @ cov_inv @ residuals

            return chi2

    # 4. 执行拟合 - 完全使用 lmfit
    success = False
    message = ""
    popt = None
    pcov = None
    lmfit_result = None

    if verbose:
        print(f"使用 lmfit 进行拟合，方法: {lmfit_method}")

    try:
        # 创建Minimizer对象
        minimizer = lmfit.Minimizer(objective, params)

        # 根据方法设置不同的参数
        fit_kws = {}

        if lmfit_method == 'least_squares':
            # 传递trf或dogbox方法给scipy的least_squares
            if method in ['trf', 'dogbox']:
                fit_kws['method'] = method
            fit_kws['ftol'] = 1e-10
            fit_kws['xtol'] = 1e-10
            fit_kws['gtol'] = 1e-10
            fit_kws['max_nfev'] = 10000 * n_params
        elif lmfit_method == 'leastsq':
            fit_kws['maxfev'] = 10000 * n_params
        elif lmfit_method == 'lbfgsb':
            # 修正：lbfgsb方法使用'tol'而不是'ftol'和'gtol'
            fit_kws['maxiter'] = 10000 * n_params
            fit_kws['tol'] = 1e-12  # 使用tol参数
        elif lmfit_method == 'differential_evolution':
            fit_kws['max_nfev'] = 10000 * n_params
            fit_kws['popsize'] = min(20, 5 * n_params)  # 种群大小
            fit_kws['tol'] = 1e-6

            # 对于differential_evolution方法，需要有限边界
            # 将无限边界替换为有限值
            for name in param_names:
                param = params[name]
                if param.min == -np.inf:
                    param.min = -1e10
                if param.max == np.inf:
                    param.max = 1e10

        elif lmfit_method == 'nelder':
            fit_kws['maxiter'] = 10000 * n_params

        # 执行最小化
        result = minimizer.minimize(method=lmfit_method, **fit_kws)
        lmfit_result = result

        if result.success:
            success = True
            message = result.message

            # 提取最优参数
            popt = np.array([result.params[name].value for name in param_names])

            if verbose:
                print(f"lmfit 拟合成功，参数: {popt}")
                print(f"消息: {message}")
        else:
            # 如果第一次失败，尝试不同的方法
            if verbose:
                print(f"初始方法 {lmfit_method} 失败，尝试备选方法...")

            # 尝试备选方法
            backup_methods = ['least_squares', 'lbfgsb', 'nelder']
            for backup_method in backup_methods:
                if backup_method != lmfit_method:
                    try:
                        if verbose:
                            print(f"尝试方法: {backup_method}")

                        # 为备选方法创建新的目标函数
                        if backup_method in ['least_squares', 'leastsq']:
                            # 需要残差向量的方法
                            def backup_objective(params):
                                param_values = np.array([params[name].value for name in param_names])
                                model_vals = model_func(x_values, param_values)
                                residuals = Y_mean - model_vals
                                try:
                                    L = cholesky(mean_cov, lower=True, check_finite=False)
                                    weighted_residuals = solve_triangular(L, residuals, lower=True)
                                except:
                                    weighted_residuals = residuals / Y_err_simple
                                return weighted_residuals
                        else:
                            # 需要标量值的方法
                            def backup_objective(params):
                                param_values = np.array([params[name].value for name in param_names])
                                model_vals = model_func(x_values, param_values)
                                residuals = Y_mean - model_vals
                                chi2 = residuals.T @ cov_inv @ residuals
                                return chi2

                        backup_minimizer = lmfit.Minimizer(backup_objective, params)

                        # 为备选方法设置参数
                        backup_fit_kws = {}
                        if backup_method == 'least_squares':
                            backup_fit_kws['ftol'] = 1e-10
                            backup_fit_kws['xtol'] = 1e-10
                            backup_fit_kws['gtol'] = 1e-10
                            backup_fit_kws['max_nfev'] = 10000 * n_params
                        elif backup_method == 'lbfgsb':
                            backup_fit_kws['maxiter'] = 10000 * n_params
                            backup_fit_kws['tol'] = 1e-12
                        elif backup_method == 'nelder':
                            backup_fit_kws['maxiter'] = 10000 * n_params

                        result = backup_minimizer.minimize(method=backup_method, **backup_fit_kws)

                        if result.success:
                            success = True
                            message = f"使用备选方法 {backup_method} 成功: {result.message}"
                            popt = np.array([result.params[name].value for name in param_names])
                            lmfit_result = result
                            break
                    except:
                        continue

            if not success:
                raise RuntimeError(
                    f"所有 lmfit 方法都失败，最后错误: {result.message if 'result' in locals() else '未知错误'}")

    except Exception as e:
        message = f"lmfit 拟合失败: {e}"
        if verbose:
            print(message)
        return {
            'params': None,
            'params_err': None,
            'params_bounds': None,
            'bounds_type': bounds_type,
            'chi2': None,
            'chi2_dof': None,
            'cov_matrix': None,
            'corr_matrix': None,
            'data_cov': mean_cov,
            'success': False,
            'message': message,
            'n_configs': n_configs,
            'n_points': m_points,
            'n_params': n_params,
            'degrees_of_freedom': None,
            'has_bounds': has_bounds,
            'residuals': None
        }

    # 5. 根据 uncertainty_method 计算参数协方差矩阵
    if success and popt is not None:
        if uncertainty_method == 'covariance':
            # 使用与 curve_fit 相同的方法计算协方差矩阵
            # 计算雅可比矩阵
            eps = 1e-8
            jacobian = np.zeros((m_points, n_params))

            for i in range(n_params):
                params_plus = popt.copy()
                params_minus = popt.copy()
                params_plus[i] += eps
                params_minus[i] -= eps

                f_plus = model_func(x_values, params_plus)
                f_minus = model_func(x_values, params_minus)
                jacobian[:, i] = (f_plus - f_minus) / (2 * eps)

            # 参数协方差 = (J^T W J)^{-1}，其中 W = C^{-1}
            try:
                Jt_W = jacobian.T @ cov_inv
                pcov = inv(Jt_W @ jacobian)
                if verbose:
                    print("使用协方差矩阵方法计算参数不确定性")
            except:
                # 如果矩阵不可逆，使用伪逆
                warnings.warn("无法计算参数协方差矩阵的逆，使用伪逆")
                Jt_W = jacobian.T @ cov_inv
                pcov = pinv(Jt_W @ jacobian)

        elif uncertainty_method == 'bootstrap':
            if verbose:
                print("使用 Bootstrap 方法估计不确定性...")

            # 简单的 Bootstrap 实现
            n_bootstrap = min(100, n_configs * 10)  # Bootstrap 样本数
            bootstrap_params = []

            for b in range(n_bootstrap):
                try:
                    # 重采样组态（有放回）
                    indices = np.random.choice(n_configs, n_configs, replace=True)
                    bootstrap_data = data_matrix[indices, :]

                    # 计算平均值和协方差
                    bootstrap_mean = np.mean(bootstrap_data, axis=0)
                    if n_configs > 1:
                        bootstrap_cov = np.cov(bootstrap_data, rowvar=False, ddof=1) / n_configs
                    else:
                        continue

                    # 对每个bootstrap样本执行相同的lmfit拟合过程
                    # 这里简化处理，使用当前最优值作为初始值
                    bootstrap_params_obj = lmfit.Parameters()
                    for i in range(n_params):
                        name = param_names[i]
                        init_val = popt[i]  # 使用当前最优值作为初始值
                        min_val, max_val = processed_bounds[i]
                        bootstrap_params_obj.add(name, value=init_val, min=min_val, max=max_val)

                    # 为bootstrap样本定义目标函数
                    def bootstrap_objective(params):
                        param_values = np.array([params[name].value for name in param_names])
                        model_vals = model_func(x_values, param_values)
                        residuals = bootstrap_mean - model_vals

                        # 使用bootstrap协方差矩阵
                        try:
                            # 计算bootstrap协方差矩阵的逆
                            bootstrap_cov_inv, _ = safe_inverse(bootstrap_cov)
                            chi2 = residuals.T @ bootstrap_cov_inv @ residuals
                            return chi2
                        except:
                            return np.sum(residuals ** 2)

                    # 最小化bootstrap样本的目标函数
                    bootstrap_minimizer = lmfit.Minimizer(bootstrap_objective, bootstrap_params_obj)

                    # 对于bootstrap，使用lbfgsb方法，因为它通常更稳定
                    bootstrap_result = bootstrap_minimizer.minimize(method='lbfgsb', maxiter=1000 * n_params, tol=1e-12)

                    if bootstrap_result.success:
                        bootstrap_popt = np.array([bootstrap_result.params[name].value for name in param_names])
                        bootstrap_params.append(bootstrap_popt)
                except:
                    # 如果拟合失败，跳过这个样本
                    continue

            if len(bootstrap_params) > 1:
                bootstrap_params = np.array(bootstrap_params)
                # 计算 Bootstrap 估计的参数协方差
                pcov = np.cov(bootstrap_params, rowvar=False, ddof=1)
                if verbose:
                    print(f"Bootstrap 完成，成功样本数: {len(bootstrap_params)}/{n_bootstrap}")
            else:
                warnings.warn("Bootstrap 方法失败，回退到协方差方法")
                # 回退到协方差方法
                eps = 1e-8
                jacobian = np.zeros((m_points, n_params))

                for i in range(n_params):
                    params_plus = popt.copy()
                    params_minus = popt.copy()
                    params_plus[i] += eps
                    params_minus[i] -= eps

                    f_plus = model_func(x_values, params_plus)
                    f_minus = model_func(x_values, params_minus)
                    jacobian[:, i] = (f_plus - f_minus) / (2 * eps)

                Jt_W = jacobian.T @ cov_inv
                pcov = pinv(Jt_W @ jacobian)

        elif uncertainty_method == 'mcmc':
            warnings.warn("MCMC 方法需要 emcee 包，暂未实现，使用协方差方法替代")
            # 使用协方差方法
            eps = 1e-8
            jacobian = np.zeros((m_points, n_params))

            for i in range(n_params):
                params_plus = popt.copy()
                params_minus = popt.copy()
                params_plus[i] += eps
                params_minus[i] -= eps

                f_plus = model_func(x_values, params_plus)
                f_minus = model_func(x_values, params_minus)
                jacobian[:, i] = (f_plus - f_minus) / (2 * eps)

            Jt_W = jacobian.T @ cov_inv
            pcov = pinv(Jt_W @ jacobian)

        else:
            raise ValueError(f"不支持的 uncertainty_method: {uncertainty_method}")

    # 6. 检查参数是否在边界上
    params_on_boundary = None
    if bounds_info is not None and popt is not None:
        params_on_boundary = []
        tolerance = 1e-6
        for i, (param, (lb, ub)) in enumerate(zip(popt, bounds_info)):
            on_lower = (lb != -np.inf) and (abs(param - lb) < tolerance)
            on_upper = (ub != np.inf) and (abs(param - ub) < tolerance)
            params_on_boundary.append(on_lower or on_upper)

            if verbose and (on_lower or on_upper):
                bound_type = "下界" if on_lower else "上界"
                bound_val = lb if on_lower else ub
                print(f"警告: 参数{i}在边界上 ({bound_type}: {bound_val:.6f})")

    # 7. 计算拟合统计量
    if popt is not None:
        residuals = Y_mean - model_func(x_values, popt)
        chi2 = residuals.T @ cov_inv @ residuals
        degrees_of_freedom = m_points - n_params

        if degrees_of_freedom <= 0:
            warnings.warn(f"自由度({degrees_of_freedom}) <= 0，χ²/dof不可靠")
            chi2_dof = float('inf') if degrees_of_freedom == 0 else float('nan')
        else:
            chi2_dof = chi2 / degrees_of_freedom

        # 计算简单χ²（忽略相关性，用于比较）
        simple_chi2 = np.sum((residuals / Y_err_simple) ** 2)
        simple_chi2_dof = simple_chi2 / degrees_of_freedom if degrees_of_freedom > 0 else np.nan
    else:
        residuals = None
        chi2 = None
        degrees_of_freedom = None
        chi2_dof = None
        simple_chi2 = None
        simple_chi2_dof = None

    # 8. 计算参数误差 - 修复警告问题
    params_err = None
    if popt is not None and pcov is not None:
        try:
            # 检查协方差矩阵对角线是否为非负数
            diag = np.diag(pcov)
            if np.all(diag >= 0) and np.all(np.isfinite(diag)):
                params_err = np.sqrt(diag)
            else:
                # 如果对角线有负数，使用绝对值或设置为NaN
                warnings.warn("协方差矩阵对角线包含负数或非有限值")
                params_err = np.full_like(popt, np.nan)
        except:
            params_err = np.full_like(popt, np.nan)
    elif popt is not None:
        params_err = np.full_like(popt, np.nan)

    # 9. 计算参数相关系数矩阵（可选）
    corr_matrix = None
    if compute_correlation and success and pcov is not None and popt is not None:
        try:
            # 相关系数矩阵 = D^{-1} * Σ * D^{-1}，其中D是对角矩阵，元素为标准差
            std_params = np.sqrt(np.diag(pcov))
            if np.all(std_params > 0) and np.all(np.isfinite(std_params)):
                D_inv = np.diag(1.0 / std_params)
                corr_matrix = D_inv @ pcov @ D_inv
                # 确保相关系数矩阵对角线为1
                np.fill_diagonal(corr_matrix, 1.0)
            else:
                corr_matrix = np.full((n_params, n_params), np.nan)
        except:
            corr_matrix = np.full((n_params, n_params), np.nan)

    # 10. 输出诊断信息
    if verbose and popt is not None:
        print("\n=== 拟合诊断信息 ===")
        print(f"自由度 (d.o.f.): {degrees_of_freedom}")
        if residuals is not None:
            print(f"残差绝对值均值: {np.mean(np.abs(residuals)):.6e}")
            print(f"残差绝对值最大值: {np.max(np.abs(residuals)):.6e}")
        print(f"χ² (关联拟合): {chi2:.6e}")
        print(f"χ²/d.o.f. (关联拟合): {chi2_dof:.6f}")
        print(f"χ² (简单，忽略相关性): {simple_chi2:.6e}")
        print(f"χ²/d.o.f. (简单): {simple_chi2_dof:.6f}")

        if degrees_of_freedom > 0 and simple_chi2_dof > 0:
            ratio = chi2_dof / simple_chi2_dof
            print(f"关联χ²/简单χ²比例: {ratio:.6f}")

            # 判断拟合质量
            if 0.5 < chi2_dof < 2.0:
                print("✅ 关联拟合质量: 良好 (0.5 < χ²/d.o.f. < 2.0)")
            elif chi2_dof < 0.5:
                print("⚠️ 关联拟合质量: 可能误差被高估或模型过拟合 (χ²/d.o.f. < 0.5)")
            else:
                print("⚠️ 关联拟合质量: 可能误差被低估或模型不合适 (χ²/d.o.f. > 2.0)")

            if ratio < 0.1:
                print("⚠️ 警告: 关联χ²远小于简单χ²，可能协方差矩阵被高估")
            elif ratio > 10:
                print("⚠️ 警告: 关联χ²远大于简单χ²，可能协方差矩阵被低估")

    # 11. 返回结果
    result_dict = {
        'params': popt,
        'params_err': params_err,
        'params_bounds': params_on_boundary,
        'bounds_type': bounds_type,
        'bounds_info': bounds_info,
        'chi2': chi2,
        'chi2_dof': chi2_dof,
        'simple_chi2': simple_chi2,
        'simple_chi2_dof': simple_chi2_dof,
        'cov_matrix': pcov,
        'corr_matrix': corr_matrix,
        'data_cov': mean_cov,
        'Y_mean': Y_mean,
        'Y_err_simple': Y_err_simple,
        'success': success,
        'message': message,
        'n_configs': n_configs,
        'n_points': m_points,
        'n_params': n_params,
        'degrees_of_freedom': degrees_of_freedom,
        'has_bounds': has_bounds,
        'residuals': residuals,
        'inv_method': inv_method,
        'uncertainty_method': uncertainty_method
    }

    # 添加 lmfit 结果对象
    if lmfit_result is not None:
        result_dict['lmfit_result'] = lmfit_result

    return result_dict