import numpy as np
import scipy.optimize as opt
from scipy.linalg import svd
import matplotlib.pyplot as plt


class JointCorrelatedFit:
    def __init__(self, data_matrices, model_functions):
        self.data_matrices = data_matrices
        self.model_functions = model_functions
        self.n_datasets = len(data_matrices)

        if len(data_matrices) != len(model_functions):
            raise ValueError("Number of data matrices must match number of model functions")

        self.calculate_joint_statistics()

    def calculate_joint_statistics(self):
        self.y_means = []
        self.cov_matrices = []
        self.inv_covs = []
        self.N_samples_list = []
        self.N_observables_list = []

        for i, data_matrix in enumerate(self.data_matrices):
            N_samples, N_observables = data_matrix.shape
            if N_samples <= 1:
                raise ValueError("At least two samples are required for a correlated fit")
            self.N_samples_list.append(N_samples)
            self.N_observables_list.append(N_observables)

            y_mean = np.mean(data_matrix, axis=0)
            cov_matrix = np.cov(data_matrix, rowvar=False, ddof=1) / N_samples
            inv_cov = self.regularized_inverse_covariance(cov_matrix)

            self.y_means.append(y_mean)
            self.cov_matrices.append(cov_matrix)
            self.inv_covs.append(inv_cov)

    def regularized_inverse_covariance(self, cov_matrix, regularization=1e-6):
        n_obs = cov_matrix.shape[0]
        regularized_cov = cov_matrix + regularization * np.eye(n_obs) * np.trace(cov_matrix) / n_obs

        try:
            U, s, Vt = svd(regularized_cov, full_matrices=False)
            max_cond = 1e10
            s_max = np.max(s)
            s_min = s_max / max_cond
            s_inv = np.array([1 / s_val if s_val > s_min else 0 for s_val in s])
            inv_cov = Vt.T @ np.diag(s_inv) @ U.T
            return inv_cov
        except:
            return np.linalg.pinv(regularized_cov)

    def joint_chi2(self, params, x_data_list=None):
        if x_data_list is None:
            x_data_list = [np.arange(N_obs) for N_obs in self.N_observables_list]

        total_chi2 = 0.0

        for i in range(self.n_datasets):
            try:
                y_data = self.y_means[i]
                inv_cov = self.inv_covs[i]
                x_data = x_data_list[i]
                model_func = self.model_functions[i]

                y_model = model_func(params, x_data)
                residuals = y_data - y_model
                chi2_i = residuals @ inv_cov @ residuals
                total_chi2 += chi2_i

            except Exception as e:
                return 1e10

        param_penalty = 1e-10 * np.sum(np.array(params) ** 2)
        total_chi2 += param_penalty

        return total_chi2

    def get_degrees_of_freedom(self):
        """Calculate total degrees of freedom"""
        total_data_points = sum(self.N_observables_list)
        n_params = len(self.fit_result.x) if hasattr(self, 'fit_result') else 0
        dof = total_data_points - n_params
        return dof, total_data_points, n_params

    def fit(self, initial_guess, x_data_list=None, bounds=None, method=None):
        if x_data_list is None:
            x_data_list = [np.arange(N_obs) for N_obs in self.N_observables_list]

        def objective(params):
            return self.joint_chi2(params, x_data_list=x_data_list)

        if method is None or len(str(method)) < 1:
            methods = ['L-BFGS-B', 'Nelder-Mead', 'Powell']
        else:
            methods = [method]

        best_result = None
        best_chi2 = np.inf

        for method in methods:
            try:
                if bounds is not None and method in ['L-BFGS-B', 'TNC', 'SLSQP']:
                    result = opt.minimize(objective, initial_guess,
                                          method=method, bounds=bounds,
                                          options={'gtol': 1e-8, 'ftol': 1e-10})
                else:
                    result = opt.minimize(objective, initial_guess,
                                          method=method,
                                          options={'gtol': 1e-8, 'ftol': 1e-10})

                if result.success and result.fun < best_chi2:
                    best_result = result
                    best_chi2 = result.fun

            except Exception as e:
                continue

        if best_result is None:
            result = opt.minimize(objective, initial_guess, method='Nelder-Mead')
            best_result = result

        self.fit_result = best_result

        # Compute parameter errors
        param_errors = self.get_parameter_errors_simple(
            x_data_list=x_data_list, best_params=best_result.x
        )
        best_result.param_errors = param_errors

        # Compute chi2 for each individual dataset
        dataset_chi2 = []
        for i in range(self.n_datasets):
            chi2_i = self._dataset_chi2(best_result.x, i, x_data_list[i] if x_data_list else None)
            dataset_chi2.append(chi2_i)
        best_result.dataset_chi2 = dataset_chi2

        # Calculate degrees of freedom information
        dof_info = self.get_degrees_of_freedom()
        best_result.dof = dof_info[0]
        best_result.total_data_points = dof_info[1]
        best_result.n_params = dof_info[2]
        best_result.chi2_per_dof = best_result.fun / best_result.dof

        return best_result

    def _dataset_chi2(self, params, dataset_idx, x_data=None):
        if x_data is None:
            x_data = np.arange(self.N_observables_list[dataset_idx])

        y_data = self.y_means[dataset_idx]
        inv_cov = self.inv_covs[dataset_idx]
        model_func = self.model_functions[dataset_idx]

        y_model = model_func(params, x_data)
        residuals = y_data - y_model
        return residuals @ inv_cov @ residuals

    def get_parameter_errors_simple(self, x_data_list=None, best_params=None):
        if best_params is None:
            if not hasattr(self, 'fit_result'):
                raise ValueError("Need to perform fit first")
            best_params = self.fit_result.x

        if x_data_list is None:
            x_data_list = [np.arange(N_obs) for N_obs in self.N_observables_list]

        n_params = len(best_params)
        best_chi2 = self.joint_chi2(best_params, x_data_list=x_data_list)

        errors = np.zeros(n_params)

        for i in range(n_params):
            def chi2_variation(value):
                test_params = best_params.copy()
                test_params[i] = value
                return self.joint_chi2(test_params, x_data_list=x_data_list) - best_chi2 - 1.0

            try:
                root_right = opt.root_scalar(chi2_variation,
                                             x0=best_params[i],
                                             x1=best_params[i] + 0.1 * abs(best_params[i]))
                if root_right.converged:
                    error_right = root_right.root - best_params[i]
                else:
                    error_right = np.nan

                root_left = opt.root_scalar(chi2_variation,
                                            x0=best_params[i],
                                            x1=best_params[i] - 0.1 * abs(best_params[i]))
                if root_left.converged:
                    error_left = best_params[i] - root_left.root
                else:
                    error_left = np.nan

                if not np.isnan(error_right) and not np.isnan(error_left):
                    errors[i] = (error_right + error_left) / 2
                elif not np.isnan(error_right):
                    errors[i] = error_right
                elif not np.isnan(error_left):
                    errors[i] = error_left
                else:
                    errors[i] = np.nan

            except:
                errors[i] = np.nan

        return errors

    def bootstrap_joint_fit(self, initial_guess, n_bootstrap=100, x_data_list=None):
        if x_data_list is None:
            x_data_list = [np.arange(N_obs) for N_obs in self.N_observables_list]

        bootstrap_params = []

        for i in range(n_bootstrap):
            try:
                bootstrap_data_matrices = []
                bootstrap_inv_covs = []
                bootstrap_y_means = []

                for j, data_matrix in enumerate(self.data_matrices):
                    n_samples = self.N_samples_list[j]
                    indices = np.random.choice(n_samples, n_samples, replace=True)
                    bootstrap_data = data_matrix[indices]

                    y_bootstrap = np.mean(bootstrap_data, axis=0)
                    cov_bootstrap = np.cov(bootstrap_data, rowvar=False, ddof=1) / n_samples
                    inv_cov_bootstrap = self.regularized_inverse_covariance(cov_bootstrap)

                    bootstrap_data_matrices.append(bootstrap_data)
                    bootstrap_y_means.append(y_bootstrap)
                    bootstrap_inv_covs.append(inv_cov_bootstrap)

                def bootstrap_chi2(params):
                    total_chi2 = 0.0
                    for k in range(self.n_datasets):
                        y_model = self.model_functions[k](params, x_data_list[k])
                        residuals = bootstrap_y_means[k] - y_model
                        total_chi2 += residuals @ bootstrap_inv_covs[k] @ residuals
                    return total_chi2

                result = opt.minimize(bootstrap_chi2, initial_guess, method='Nelder-Mead')
                if result.success:
                    bootstrap_params.append(result.x)

            except Exception as e:
                continue

        if len(bootstrap_params) == 0:
            print("Warning: All bootstrap iterations failed")
            return np.full(len(initial_guess), np.nan), np.array([])

        bootstrap_params = np.array(bootstrap_params)
        bootstrap_errors = np.std(bootstrap_params, axis=0)

        return bootstrap_errors, bootstrap_params

    def print_joint_fit_results(self, param_names=None, true_values=None):
        """
        Print formatted joint fit results with proper statistical information
        """
        if not hasattr(self, 'fit_result'):
            raise ValueError("Need to perform fit first")

        result = self.fit_result
        n_params = len(result.x)

        if param_names is None:
            param_names = [f'Param_{i}' for i in range(n_params)]

        # Get degrees of freedom information
        dof = result.dof
        total_data_points = result.total_data_points
        n_params = result.n_params
        chi2_per_dof = result.chi2_per_dof

        print("\n" + "=" * 60)
        print("JOINT CORRELATED FIT RESULTS")
        print("=" * 60)

        print(f"Success: {result.success}")
        print(f"Message: {result.message}")

        print(f"\nFit Quality Statistics:")
        print(f"  Total χ²: {result.fun:.4f}")
        print(f"  Degrees of freedom: {dof}")
        print(f"  χ²/dof: {chi2_per_dof:.4f}")
        print(f"  Data points: {total_data_points}")
        print(f"  Parameters: {n_params}")

        # Interpret χ²/dof
        if chi2_per_dof < 0.5:
            quality = "EXCELLENT (might be overfitting)"
        elif 0.5 <= chi2_per_dof <= 1.5:
            quality = "GOOD"
        elif 1.5 < chi2_per_dof <= 2.0:
            quality = "ACCEPTABLE"
        elif 2.0 < chi2_per_dof <= 3.0:
            quality = "POOR"
        else:
            quality = "BAD (model might be wrong or errors underestimated)"

        print(f"  Quality: {quality}")

        # Print chi2 for each dataset
        if hasattr(result, 'dataset_chi2'):
            print(f"\nχ² per dataset:")
            for i, chi2 in enumerate(result.dataset_chi2):
                dataset_dof = self.N_observables_list[i]
                chi2_dof_i = chi2 / dataset_dof
                print(f"  Dataset {i}: χ² = {chi2:.4f}, dof = {dataset_dof}, χ²/dof = {chi2_dof_i:.4f}")

        print(f"\nFitted parameters with errors:")
        print("-" * 70)

        if hasattr(result, 'param_errors'):
            errors = result.param_errors
            for i, name in enumerate(param_names):
                line = f"{name:10s} = {result.x[i]:10.6f} ± {errors[i]:8.6f}"
                if true_values is not None and i < len(true_values):
                    true_val = true_values[i]
                    diff = abs(result.x[i] - true_val)
                    diff_sigma = diff / errors[i] if errors[i] > 0 else np.inf
                    line += f"  (true: {true_val:.6f}, diff: {diff:.6f}, {diff_sigma:.1f}σ)"
                print(line)
        else:
            for i, name in enumerate(param_names):
                line = f"{name:10s} = {result.x[i]:10.6f}"
                if true_values is not None and i < len(true_values):
                    true_val = true_values[i]
                    diff = abs(result.x[i] - true_val)
                    line += f"  (true: {true_val:.6f}, diff: {diff:.6f})"
                print(line)

    def plot_joint_fit(self, x_data_list=None, dataset_names=None, true_params=None):
        if not hasattr(self, 'fit_result'):
            raise ValueError("Need to perform fit first")

        if x_data_list is None:
            x_data_list = [np.arange(N_obs) for N_obs in self.N_observables_list]

        if dataset_names is None:
            dataset_names = [f'Dataset {i}' for i in range(self.n_datasets)]

        n_cols = min(3, self.n_datasets)
        n_rows = (self.n_datasets + n_cols - 1) // n_cols

        plt.figure(figsize=(5 * n_cols, 4 * n_rows))

        for i in range(self.n_datasets):
            plt.subplot(n_rows, n_cols, i + 1)

            y_mean = self.y_means[i]
            y_std = np.std(self.data_matrices[i], axis=0)
            x_data = x_data_list[i]

            y_fit = self.model_functions[i](self.fit_result.x, x_data)

            plt.errorbar(x_data, y_mean, yerr=y_std, fmt='o', alpha=0.7,
                         label=f'{dataset_names[i]} data')
            plt.plot(x_data, y_fit, 'r-', linewidth=2, label='Joint fit')

            if true_params is not None:
                y_true = self.model_functions[i](true_params, x_data)
                plt.plot(x_data, y_true, 'g--', label='True')

            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend()

            # Show χ²/dof for this dataset in the title
            if hasattr(self.fit_result, 'dataset_chi2'):
                chi2 = self.fit_result.dataset_chi2[i]
                dof_i = self.N_observables_list[i]
                chi2_dof_i = chi2 / dof_i
                plt.title(f'{dataset_names[i]}\nχ²/dof = {chi2_dof_i:.3f}')
            else:
                plt.title(f'{dataset_names[i]}')

        plt.tight_layout()
        plt.show()


# Example usage with improved statistics display
def demo_joint_fit():
    # Shared true parameters
    true_A, true_B, true_tau = 2.5, 1.0, 3.0
    true_params = [true_A, true_B, true_tau]
    param_names = ['A', 'B', 'tau']

    # Define different model functions
    def model1(params, x):
        A, B, tau = params
        return A * np.exp(-x / tau) + B

    def model2(params, x):
        A, B, tau = params
        return A * (1 - np.exp(-x / tau)) + B

    def model3(params, x):
        A, B, tau = params
        return A * np.exp(-(x / tau) ** 2) + B

    model_functions = [model1, model2, model3]

    # Generate data for each dataset
    np.random.seed(42)
    data_matrices = []
    x_data_list = []

    for i, model_func in enumerate(model_functions):
        N_samples = 200 + i * 50
        N_points = 15 + i * 5

        x_data = np.linspace(0, 10, N_points)
        x_data_list.append(x_data)

        noise_cov = 0.1 * np.exp(-0.2 * np.abs(np.subtract.outer(x_data, x_data)))
        correlated_noise = np.random.multivariate_normal(
            np.zeros(N_points), noise_cov, N_samples
        )

        data_matrix = np.zeros((N_samples, N_points))
        for j in range(N_samples):
            true_signal = model_func(true_params, x_data)
            data_matrix[j] = true_signal + correlated_noise[j]

        data_matrices.append(data_matrix)
        print(f"Dataset {i}: {data_matrix.shape} (samples × observables)")

    # Create joint fit object
    joint_fitter = JointCorrelatedFit(data_matrices, model_functions)

    # Initial guess
    initial_guess = [2.0, 1.2, 2.5]

    # Perform joint fit
    result = joint_fitter.fit(initial_guess, x_data_list=x_data_list)

    # Print results with improved statistics
    joint_fitter.print_joint_fit_results(param_names=param_names, true_values=true_params)

    # Plot results
    dataset_names = ['Exp Decay', 'Modified Exp', 'Gaussian Exp']
    joint_fitter.plot_joint_fit(
        x_data_list=x_data_list,
        dataset_names=dataset_names,
        true_params=true_params
    )

    return joint_fitter, result

def JointCorrelatedFitAndDraw(xdata, ydata, func, guess, bound=None, method=""):
    """

    :param xdata: [xlst1, xlst2, ...], xlst: arraylike
    :param ydata: [ymatrix1, ymatrix2, ...]
    :param func: [func1, func2, ... ]
    :param guess: arraylike
    :param bound: [[min1, max1], [min2, max2], ...]
    :param method:
    """
    # Create joint fit object
    joint_fitter = JointCorrelatedFit(ydata, func)

    # Perform joint fit
    joint_fitter.fit(guess, x_data_list=xdata, bounds=bound, method=method)

    # Print results with improved statistics
    joint_fitter.print_joint_fit_results()

    # Plot results
    joint_fitter.plot_joint_fit(x_data_list=xdata)

if __name__ == "__main__":
    joint_fitter, result = demo_joint_fit()
