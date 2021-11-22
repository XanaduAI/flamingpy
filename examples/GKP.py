import numpy as np
from ft_stack import GKP, viz

alpha = np.sqrt(np.pi)
xs = np.arange(-10, 10, 0.01)

ns, fs = GKP.integer_fractional(xs, alpha)
viz.plot_integer_fractional(xs, ns, fs, alpha)

bit_values = GKP.GKP_binner(xs)
viz.plot_GKP_bins(xs, bit_values, alpha)

delta = 0.1
error_hom_val = GKP.Z_err_cond([delta] * len(xs), xs, use_hom_val=True)    
error_no_hom_val = GKP.Z_err_cond([delta] * len(xs), xs)

viz.plt_Z_err_cond(xs, error_hom_val, alpha, True)
viz.plt_Z_err_cond(xs, error_no_hom_val, alpha, False)
