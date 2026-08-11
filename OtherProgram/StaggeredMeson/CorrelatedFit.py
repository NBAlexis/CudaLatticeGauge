import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.linalg import inv, pinv, cholesky
import warnings
from typing import Optional, Tuple, List, Callable, Union


def correlated_fit_with_flexible_bounds(
        data_matrix: np.ndarray,
        x_values: np.ndarray,
        model_func: Callable,
        initial_params: np.ndarray,
        param_bounds: Optional[
            Union[List[Tuple[Optional[float], Optional[float]]], Tuple[List[float], List[float]]]] = None,
        method: str = 'trf',
        compute_correlation: bool = True,
        verbose: bool = False
) -> dict:
    """
    执行关联拟合，支持灵活的边界约束

    关键修正：正确计算平均值的协方差矩阵

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
        优化方法：
        - 'trf': 信赖域反射法（推荐，支持边界）
        - 'dogbox': dogbox算法（支持边界）
        - 'lm': Levenberg-Marquardt（不支持边界，无约束时使用）
    compute_correlation : bool
        是否计算参数间的相关系数矩阵
    verbose : bool
        是否打印详细过程信息

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

    # 处理边界 - 支持多种格式
    bounds_for_scipy = None
    bounds_info = None
    bounds_type = None
    has_bounds = False

    if param_bounds is not None:
        # 检查边界格式
        if isinstance(param_bounds, tuple) and len(param_bounds) == 2:
            # 传统格式: ([min1, min2, ...], [max1, max2, ...])
            lower_list, upper_list = param_bounds

            if len(lower_list) != n_params or len(upper_list) != n_params:
                raise ValueError(
                    f"传统边界格式: 下界列表长度({len(lower_list)})和上界列表长度({len(upper_list)})必须与参数个数({n_params})相同")

            # 转换为灵活格式
            param_bounds_flexible = []
            for i in range(n_params):
                min_val = lower_list[i]
                max_val = upper_list[i]

                # 将np.inf/np.nan转换为None
                if min_val is None or np.isnan(min_val) or np.isneginf(min_val):
                    min_val = None
                elif np.isinf(min_val):
                    min_val = -np.inf if min_val < 0 else np.inf

                if max_val is None or np.isnan(max_val) or np.isposinf(max_val):
                    max_val = None
                elif np.isinf(max_val):
                    max_val = np.inf if max_val > 0 else -np.inf

                param_bounds_flexible.append((min_val, max_val))

            param_bounds = param_bounds_flexible

        elif isinstance(param_bounds, list):
            # 灵活格式: [(min1, max1), (min2, max2), ...]
            if len(param_bounds) != n_params:
                raise ValueError(f"灵活边界格式: 边界列表长度({len(param_bounds)})必须与参数个数({n_params})相同")

            # 确保每个元素都是元组
            param_bounds_flexible = []
            for i, bound in enumerate(param_bounds):
                if isinstance(bound, tuple) and len(bound) == 2:
                    param_bounds_flexible.append(bound)
                else:
                    raise ValueError(f"边界元素{i}不是有效的元组格式: {bound}")

            param_bounds = param_bounds_flexible
        else:
            raise ValueError(
                f"不支持的边界格式: {type(param_bounds)}。应为None、元组([min...],[max...])或列表[(min,max),...]")

        # 现在param_bounds是灵活格式的列表
        # 转换边界为SciPy需要的格式
        lower_bounds = []
        upper_bounds = []
        bounds_info = []
        bounds_type = []

        has_bounds = False
        for i, (min_val, max_val) in enumerate(param_bounds):
            # 处理min_val
            if min_val is None:
                lower_bounds.append(-np.inf)
            elif np.isinf(min_val) and min_val < 0:
                lower_bounds.append(-np.inf)
            else:
                lower_bounds.append(float(min_val))
                if not np.isinf(min_val):
                    has_bounds = True

            # 处理max_val
            if max_val is None:
                upper_bounds.append(np.inf)
            elif np.isinf(max_val) and max_val > 0:
                upper_bounds.append(np.inf)
            else:
                upper_bounds.append(float(max_val))
                if not np.isinf(max_val):
                    has_bounds = True

            # 记录边界信息
            bounds_info.append((lower_bounds[-1], upper_bounds[-1]))

            # 记录边界类型
            lb_finite = not np.isinf(lower_bounds[-1]) and lower_bounds[-1] != -np.inf
            ub_finite = not np.isinf(upper_bounds[-1]) and upper_bounds[-1] != np.inf

            if not lb_finite and not ub_finite:
                bounds_type.append("无限制")
            elif lb_finite and ub_finite:
                bounds_type.append("双边限制")
            elif lb_finite:
                bounds_type.append("下界限制")
            else:
                bounds_type.append("上界限制")

        if has_bounds:
            bounds_for_scipy = (lower_bounds, upper_bounds)

            # 对于有边界的情况，强制使用支持边界的方法
            if method == 'lm':
                warnings.warn("方法'lm'不支持边界，已自动切换到'trf'")
                method = 'trf'

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

    if verbose:
        print(f"数据统计: {n_configs}个组态, {m_points}个测量点, {n_params}个参数")
        if bounds_type is not None:
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

    # 检查协方差矩阵是否可逆
    try:
        # 尝试Cholesky分解验证正定性
        L = cholesky(mean_cov, lower=True)
        cov_inv = inv(mean_cov)
        if verbose:
            print("协方差矩阵正定，使用标准逆矩阵")
    except:
        # 如果不正定，使用伪逆并发出警告
        warnings.warn("协方差矩阵不正定，使用伪逆。结果可能不稳定。")
        cov_inv = pinv(mean_cov)
        if verbose:
            print("使用伪逆矩阵")

    # 2. 定义带权重的χ²函数
    def chi2_function(params):
        """计算χ²值"""
        model_values = model_func(x_values, params)
        residuals = Y_mean - model_values

        # 使用完整的协方差矩阵：χ² = r^T C^{-1} r
        chi2 = residuals @ cov_inv @ residuals
        return chi2

    # 3. 执行拟合（支持边界）
    success = False
    message = ""
    popt = None
    pcov = None

    # 定义包装函数
    def wrapped_model_func(x, *params):
        return model_func(x, params)

    try:
        # 方法1：使用curve_fit
        if verbose:
            print(f"尝试使用curve_fit进行拟合，方法: {method}")
            if bounds_for_scipy is not None:
                print(f"使用边界: {bounds_for_scipy}")

        # 执行curve_fit - 关键：使用mean_cov作为sigma
        if bounds_for_scipy is not None:
            popt, pcov = curve_fit(
                wrapped_model_func,
                x_values,
                Y_mean,
                p0=initial_params,
                sigma=mean_cov,  # 关键修正：使用平均值的协方差
                absolute_sigma=True,  # 关键：使用绝对协方差
                method=method,
                bounds=bounds_for_scipy,
                maxfev=10000 * n_params
            )
        else:
            popt, pcov = curve_fit(
                wrapped_model_func,
                x_values,
                Y_mean,
                p0=initial_params,
                sigma=mean_cov,  # 关键修正：使用平均值的协方差
                absolute_sigma=True,  # 关键：使用绝对协方差
                method=method,
                maxfev=10000 * n_params
            )

        success = True
        message = "curve_fit成功"

        if verbose:
            print(f"curve_fit成功，参数: {popt}")

    except Exception as e_curve:
        # 方法2：如果curve_fit失败，使用minimize方法
        if verbose:
            print(f"curve_fit失败: {e_curve}，尝试minimize方法")

        try:
            # 准备minimize的边界
            if bounds_for_scipy is not None:
                minimize_bounds = []
                for i, (lb, ub) in enumerate(zip(bounds_for_scipy[0], bounds_for_scipy[1])):
                    minimize_bounds.append((lb, ub))
            else:
                minimize_bounds = None

            # 执行最小化
            result = minimize(
                chi2_function,
                initial_params,
                method='L-BFGS-B',  # L-BFGS-B支持边界
                bounds=minimize_bounds,
                options={'maxiter': 1000, 'disp': verbose, 'ftol': 1e-9}
            )

            if not result.success:
                message = f"最小化失败: {result.message}"
                raise RuntimeError(message)

            popt = result.x
            success = True
            message = f"minimize成功: {result.message}"

            if verbose:
                print(f"minimize成功，参数: {popt}")

            # 计算参数的协方差矩阵
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

            # 关键修正：正确的参数协方差矩阵计算公式
            # 对于加权最小二乘，参数协方差 = (J^T W J)^{-1}
            # 其中 W = C^{-1} 是平均值的协方差矩阵的逆
            try:
                Jt_W = jacobian.T @ cov_inv
                param_cov = inv(Jt_W @ jacobian)
                pcov = param_cov
            except:
                # 如果矩阵不可逆，使用伪逆
                warnings.warn("无法计算参数协方差矩阵的逆，使用伪逆")
                Jt_W = jacobian.T @ cov_inv
                param_cov = pinv(Jt_W @ jacobian)
                pcov = param_cov

        except Exception as e2:
            message = f"所有拟合方法都失败: {e2}"
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
                'data_cov': mean_cov,  # 返回平均值的协方差
                'success': False,
                'message': message,
                'n_configs': n_configs,
                'n_points': m_points,
                'n_params': n_params,
                'degrees_of_freedom': None,
                'has_bounds': has_bounds,
                'residuals': None
            }

    # 4. 检查参数是否在边界上
    params_on_boundary = None
    if bounds_info is not None:
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

    # 5. 计算拟合统计量
    residuals = Y_mean - model_func(x_values, popt)

    # 关键修正：使用正确的协方差矩阵计算χ²
    chi2 = residuals @ cov_inv @ residuals

    # 自由度 = 数据点数 - 参数个数
    degrees_of_freedom = m_points - n_params

    if degrees_of_freedom <= 0:
        warnings.warn(f"自由度({degrees_of_freedom}) <= 0，χ²/dof不可靠")
        chi2_dof = float('inf') if degrees_of_freedom == 0 else float('nan')
    else:
        chi2_dof = chi2 / degrees_of_freedom

    # 6. 计算简单χ²（忽略相关性，用于比较）
    simple_chi2 = np.sum((residuals / Y_err_simple) ** 2)
    simple_chi2_dof = simple_chi2 / degrees_of_freedom if degrees_of_freedom > 0 else np.nan

    # 7. 计算参数误差
    try:
        if pcov is not None and np.all(np.isfinite(pcov)):
            params_err = np.sqrt(np.diag(pcov))
        else:
            params_err = np.full_like(popt, np.nan)
    except:
        params_err = np.full_like(popt, np.nan)

    # 8. 计算参数相关系数矩阵（可选）
    corr_matrix = None
    if compute_correlation and success and pcov is not None:
        try:
            # 相关系数矩阵 = D^{-1} * Σ * D^{-1}，其中D是对角矩阵，元素为标准差
            std_params = np.sqrt(np.diag(pcov))
            if np.all(std_params > 0) and np.all(np.isfinite(std_params)):
                D_inv = np.diag(1.0 / std_params)
                corr_matrix = D_inv @ pcov @ D_inv
            else:
                corr_matrix = np.full((n_params, n_params), np.nan)
        except:
            corr_matrix = np.full((n_params, n_params), np.nan)

    # 9. 输出诊断信息
    if verbose:
        print("\n=== 拟合诊断信息 ===")
        print(f"自由度 (d.o.f.): {degrees_of_freedom}")
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

    # 10. 返回结果
    return {
        'params': popt,
        'params_err': params_err,
        'params_bounds': params_on_boundary,
        'bounds_type': bounds_type,
        'bounds_info': bounds_info,
        'chi2': chi2,
        'chi2_dof': chi2_dof,
        'simple_chi2': simple_chi2,  # 新增：简单χ²（忽略相关性）
        'simple_chi2_dof': simple_chi2_dof,  # 新增：简单χ²/d.o.f.
        'cov_matrix': pcov,
        'corr_matrix': corr_matrix,
        'data_cov': mean_cov,  # 关键修正：返回平均值的协方差
        'Y_mean': Y_mean,  # 新增：平均值
        'Y_err_simple': Y_err_simple,  # 新增：简单标准误差
        'success': success,
        'message': message,
        'n_configs': n_configs,
        'n_points': m_points,
        'n_params': n_params,
        'degrees_of_freedom': degrees_of_freedom,
        'has_bounds': has_bounds,
        'residuals': residuals
    }


# 测试函数
def test_corrected_fit():
    """测试修正后的函数"""
    np.random.seed(42)

    # 定义测试模型
    def test_model(x, params):
        a, b, c = params
        return a * np.exp(-b * x) + c

    # 生成测试数据
    n_configs = 100
    n_points = 20
    x = np.linspace(0, 5, n_points)

    # 真实参数
    true_params = [2.0, 0.5, 0.1]

    # 生成相关数据
    data = np.zeros((n_configs, n_points))
    for i in range(n_configs):
        # 共同噪声
        common_noise = np.random.normal(0, 0.5)

        # 真实值
        y_true = test_model(x, true_params)

        # 相关噪声和独立噪声
        correlated_noise = common_noise * np.ones_like(x) * 0.3
        independent_noise = np.random.normal(0, 0.2, size=n_points)

        data[i, :] = y_true + correlated_noise + independent_noise

    # 边界
    lower_bounds = [0.0, 0.0, -1.0]
    upper_bounds = [5.0, 2.0, 2.0]

    print("=" * 70)
    print("测试修正后的关联拟合函数")
    print("=" * 70)

    # 执行拟合
    result = correlated_fit_with_flexible_bounds(
        data_matrix=data,
        x_values=x,
        model_func=test_model,
        initial_params=[1.5, 0.3, 0.0],
        param_bounds=(lower_bounds, upper_bounds),
        verbose=True
    )

    if result['success']:
        print("\n拟合结果:")
        print(f"参数值: {result['params']}")
        print(f"参数误差: {result['params_err']}")
        print(f"χ²/d.o.f. (关联): {result['chi2_dof']:.4f}")
        print(f"χ²/d.o.f. (简单): {result['simple_chi2_dof']:.4f}")

        # 检查参数是否在边界上
        if result['params_bounds'] is not None:
            on_boundary = [i for i, on_b in enumerate(result['params_bounds']) if on_b]
            if on_boundary:
                print(f"在边界上的参数: {on_boundary}")

    return result


if __name__ == "__main__":
    test_corrected_fit()