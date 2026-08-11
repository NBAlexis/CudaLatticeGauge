import numpy as np
import scipy.optimize as opt
from scipy.linalg import inv, pinv, svd
import matplotlib.pyplot as plt


class RobustCorrelatedFit:
    def __init__(self, data_matrix, model_function):
        """
        Robust Correlated Fit Class

        Parameters:
        -----------
        data_matrix : ndarray, shape (N_samples, N_observables)
            Data matrix, each row is a sample, each column is an observable
        model_function : callable
            Model function, format: f(params, x) -> y
        """
        self.data_matrix = data_matrix
        self.model_function = model_function
        self.N_samples, self.N_observables = data_matrix.shape

        # Calculate statistics with improved numerical stability
        self.calculate_robust_statistics()

    def calculate_robust_statistics(self):
        """Calculate data statistics with improved numerical stability"""
        self.y_mean = np.mean(self.data_matrix, axis=0)
        if self.N_samples <= 1:
            raise ValueError("At least two samples are required for a correlated fit")
        self.sample_cov_matrix = np.cov(self.data_matrix, rowvar=False, ddof=1)
        self.cov_matrix = self.sample_cov_matrix / self.N_samples

        # Improved covariance matrix regularization
        self.inv_cov = self.regularized_inverse_covariance()

    def regularized_inverse_covariance(self, regularization=1e-6):
        """
        Compute regularized inverse covariance matrix for better numerical stability
        """
        # Add small regularization to diagonal to improve condition number
        n_obs = self.N_observables
        regularized_cov = self.cov_matrix + regularization * np.eye(n_obs) * np.trace(self.cov_matrix) / n_obs

        # Use SVD for stable inversion
        try:
            U, s, Vt = svd(regularized_cov, full_matrices=False)
            # Condition number threshold
            max_cond = 1e10
            s_max = np.max(s)
            s_min = s_max / max_cond
            s_inv = np.array([1 / s_val if s_val > s_min else 0 for s_val in s])
            inv_cov = Vt.T @ np.diag(s_inv) @ U.T
            return inv_cov
        except:
            print("SVD failed, using pseudo-inverse")
            return pinv(regularized_cov)

    def chi2(self, params, x_data=None, y_data=None, inv_cov=None):
        """
        Calculate chi2 function considering correlations
        """
        if y_data is None:
            y_data = self.y_mean
        if inv_cov is None:
            inv_cov = self.inv_cov
        if x_data is None:
            x_data = np.arange(len(y_data))

        try:
            # Calculate model prediction
            y_model = self.model_function(params, x_data)

            # Calculate residuals
            residuals = y_data - y_model

            # Mahalanobis distance: residuals^T * inv_cov * residuals
            chi2_value = residuals @ inv_cov @ residuals

            # Add small penalty for extreme parameter values to improve stability
            param_penalty = 1e-10 * np.sum(np.array(params) ** 2)
            chi2_value += param_penalty

            return chi2_value
        except (ValueError, TypeError, RuntimeError) as e:
            # Return large chi2 for invalid parameters
            return 1e10

    def fit(self, initial_guess, x_data=None, bounds=None, method=None):
        """
        Perform robust correlated fit with multiple optimization strategies
        """
        if x_data is None:
            x_data = np.arange(self.N_observables)

        # Define objective function
        def objective(params):
            return self.chi2(params, x_data=x_data)

        # If no method specified, try multiple methods
        if method is None or len(str(method)) < 1:
            methods = ['L-BFGS-B', 'Nelder-Mead', 'Powell']
        else:
            methods = [method]

        best_result = None
        best_chi2 = np.inf

        for method in methods:
            try:
                if bounds is not None and method in ['L-BFGS-B', 'TNC', 'SLSQP']:
                    print(f'fit using : {method}')
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
                print(f"Method {method} failed: {e}")
                continue

        if best_result is None:
            # If all methods fail, use the first one that didn't crash
            result = opt.minimize(objective, initial_guess, method='Nelder-Mead')
            best_result = result

        self.fit_result = best_result
        return best_result

    def bootstrap_fit(self, initial_guess, n_bootstrap=100, x_data=None):
        """
        Improved Bootstrap method with better error handling
        """
        if x_data is None:
            x_data = np.arange(self.N_observables)

        bootstrap_params = []
        n_samples = self.N_samples

        for i in range(n_bootstrap):
            # Resample with replacement
            indices = np.random.choice(n_samples, n_samples, replace=True)
            bootstrap_data = self.data_matrix[indices]

            try:
                # Calculate statistics for bootstrap sample
                y_bootstrap = np.mean(bootstrap_data, axis=0)
                cov_bootstrap = np.cov(bootstrap_data, rowvar=False, ddof=1) / n_samples

                # Regularize the bootstrap covariance matrix
                n_obs = self.N_observables
                regularized_cov = cov_bootstrap + 1e-6 * np.eye(n_obs) * np.trace(cov_bootstrap) / n_obs

                try:
                    inv_cov_bootstrap = inv(regularized_cov)
                except np.linalg.LinAlgError:
                    inv_cov_bootstrap = pinv(regularized_cov)

                # Fit bootstrap sample
                def bootstrap_chi2(params):
                    y_model = self.model_function(params, x_data)
                    residuals = y_bootstrap - y_model
                    return residuals @ inv_cov_bootstrap @ residuals

                result = opt.minimize(bootstrap_chi2, initial_guess, method='Nelder-Mead')
                if result.success:
                    bootstrap_params.append(result.x)

            except Exception as e:
                # Skip failed bootstrap iterations
                continue

        if len(bootstrap_params) == 0:
            print("Warning: All bootstrap iterations failed")
            return np.full(len(initial_guess), np.nan), np.array([])

        bootstrap_params = np.array(bootstrap_params)
        bootstrap_errors = np.std(bootstrap_params, axis=0)

        return bootstrap_errors, bootstrap_params

    def residual_analysis(self, x_data=None):
        """
        Analyze residuals to diagnose fit quality
        """
        if not hasattr(self, 'fit_result'):
            raise ValueError("Need to perform fit first")

        if x_data is None:
            x_data = np.arange(self.N_observables)

        best_params = self.fit_result.x
        y_mean = self.y_mean
        y_model = self.model_function(best_params, x_data)
        residuals = y_mean - y_model

        # Calculate normalized residuals using covariance matrix
        try:
            # Cholesky decomposition for Mahalanobis transformation
            L = np.linalg.cholesky(self.cov_matrix)
            normalized_residuals = np.linalg.solve(L, residuals)
        except:
            # Fallback: use diagonal elements only
            variances = np.diag(self.cov_matrix)
            normalized_residuals = residuals / np.sqrt(variances)

        return {
            'residuals': residuals,
            'normalized_residuals': normalized_residuals,
            'residual_std': np.std(residuals),
            'residual_mean': np.mean(residuals),
            'chi2_per_dof': self.fit_result.fun / (len(y_mean) - len(best_params))
        }


# Example with tanh function
def demo_tanh_fit():
    # Generate example data with tanh function
    np.random.seed(42)

    # True parameters for tanh function: f(x) = A * tanh((x - x0)/w) + B
    true_A, true_B, true_x0, true_w = 2.5, 1.0, 5.0, 1.5

    # Tanh model function
    def tanh_model(params, x):
        A, B, x0, w = params
        return A * np.tanh((x - x0) / w) + B

    # Generate data
    N_samples, N_points = 300, 20
    x_data = np.linspace(0, 10, N_points)

    # Generate correlated noise
    noise_cov = 0.1 * np.exp(-0.2 * np.abs(np.subtract.outer(x_data, x_data)))
    correlated_noise = np.random.multivariate_normal(
        np.zeros(N_points), noise_cov, N_samples
    )

    # Generate data matrix
    data_matrix = np.zeros((N_samples, N_points))
    for i in range(N_samples):
        true_signal = tanh_model([true_A, true_B, true_x0, true_w], x_data)
        data_matrix[i] = true_signal + correlated_noise[i]

    print(f"Data matrix shape: {data_matrix.shape}")

    # Create robust correlated fit object
    fitter = RobustCorrelatedFit(data_matrix, tanh_model)

    # Initial guess (close to true values as mentioned)
    initial_guess = [2.0, 1.2, 4.5, 1.2]  # Close to [2.5, 1.0, 5.0, 1.5]

    # Perform fit
    result = fitter.fit(initial_guess, x_data=x_data)

    print("\n=== Fit Results ===")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Fitted parameters: {result.x}")
    print(f"True parameters: [{true_A}, {true_B}, {true_x0}, {true_w}]")
    print(f"Chi2: {result.fun:.4f}")

    # Residual analysis
    residuals_info = fitter.residual_analysis(x_data=x_data)
    print(f"\nResidual Analysis:")
    print(f"Residual mean: {residuals_info['residual_mean']:.4f}")
    print(f"Residual std: {residuals_info['residual_std']:.4f}")
    print(f"Chi2/dof: {residuals_info['chi2_per_dof']:.4f}")

    # Bootstrap error analysis
    bootstrap_errors, bootstrap_params = fitter.bootstrap_fit(
        initial_guess, n_bootstrap=100, x_data=x_data
    )
    print(f"\nParameter errors (Bootstrap): {bootstrap_errors}")

    # Plot results
    plt.figure(figsize=(12, 5))

    plt.subplot(131)
    y_mean = np.mean(data_matrix, axis=0)
    y_std = np.std(data_matrix, axis=0)

    plt.errorbar(x_data, y_mean, yerr=y_std, fmt='o', alpha=0.7, label='Data')
    y_fit = tanh_model(result.x, x_data)
    plt.plot(x_data, y_fit, 'r-', linewidth=2, label='Fit')
    plt.plot(x_data, tanh_model([true_A, true_B, true_x0, true_w], x_data),
             'g--', label='True')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Tanh Function Fit')

    plt.subplot(132)
    residuals = residuals_info['residuals']
    plt.errorbar(x_data, residuals, yerr=y_std, fmt='o')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('Residuals')
    plt.title('Residuals')

    plt.subplot(133)
    param_names = ['A', 'B', 'x0', 'w']
    for i, param_name in enumerate(param_names):
        plt.hist(bootstrap_params[:, i], alpha=0.7, label=param_name, bins=15)
    plt.axvline(x=true_A, color='r', linestyle='--', label='True values')
    plt.axvline(x=true_B, color='r', linestyle='--')
    plt.axvline(x=true_x0, color='r', linestyle='--')
    plt.axvline(x=true_w, color='r', linestyle='--')
    plt.xlabel('Parameter value')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title('Bootstrap Parameter Distribution')

    plt.tight_layout()
    plt.show()

    # Additional diagnostic plots
    plt.figure(figsize=(10, 4))

    plt.subplot(121)
    # Normalized residuals Q-Q plot
    normalized_residuals = residuals_info['normalized_residuals']
    plt.scatter(np.sort(normalized_residuals),
                np.sort(np.random.normal(0, 1, len(normalized_residuals))),
                alpha=0.7)
    plt.plot([-3, 3], [-3, 3], 'r--')
    plt.xlabel('Normalized Residuals')
    plt.ylabel('Theoretical Normal Quantiles')
    plt.title('Q-Q Plot for Residuals')

    plt.subplot(122)
    # Correlation matrix
    corr_matrix = fitter.cov_matrix / np.outer(np.sqrt(np.diag(fitter.cov_matrix)),
                                               np.sqrt(np.diag(fitter.cov_matrix)))
    plt.imshow(corr_matrix, cmap='coolwarm', aspect='auto', vmin=-1, vmax=1)
    plt.colorbar(label='Correlation')
    plt.title('Data Correlation Matrix')

    plt.tight_layout()
    plt.show()

    return fitter, result

def CorrelatedFitAndDraw(xdata, ydata, func, guess, bound=None, method=""):
    # Create robust correlated fit object
    fitter = RobustCorrelatedFit(ydata, func)
    # Perform fit
    if len(method) < 1 or bound is None:
        result = fitter.fit(guess, x_data=xdata)
    else:
        result = fitter.fit(guess, x_data=xdata, bounds=bound, method=method)
    print("\n=== Fit Results ===")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    # print(f"Chi2: {result.fun:.4f}")

    # Residual analysis
    residuals_info = fitter.residual_analysis(x_data=xdata)
    # print(f"\nResidual Analysis:")
    # print(f"Residual mean: {residuals_info['residual_mean']:.4f}")
    # print(f"Residual std: {residuals_info['residual_std']:.4f}")
    print(f"Chi2/dof: {residuals_info['chi2_per_dof']:.4f}")

    # Bootstrap error analysis
    bootstrap_errors, bootstrap_params = fitter.bootstrap_fit(
        guess, n_bootstrap=100, x_data=xdata
    )
    print(f"Fitted parameters: {result.x}, Parameter errors (Bootstrap): {bootstrap_errors}")

    # Plot results
    plt.figure(figsize=(12, 5))

    plt.subplot(131)
    y_mean = np.mean(ydata, axis=0)
    y_std = np.std(ydata, axis=0)

    plt.errorbar(xdata, y_mean, yerr=y_std, fmt='o', alpha=0.7, label='Data')
    y_fit = func(result.x, xdata)
    plt.plot(xdata, y_fit, 'r-', linewidth=2, label='Fit')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Tanh Function Fit')

    plt.subplot(132)
    residuals = residuals_info['residuals']
    plt.errorbar(xdata, residuals, yerr=y_std, fmt='o')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('x')
    plt.ylabel('Residuals')
    plt.title('Residuals')

    plt.subplot(133)
    for i in range(len(guess)):
        plt.hist(bootstrap_params[:, i], alpha=0.7, bins=15)
    plt.xlabel('Parameter value')
    plt.ylabel('Frequency')
    plt.title('Bootstrap Parameter Distribution')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    fitter, result = demo_tanh_fit()
