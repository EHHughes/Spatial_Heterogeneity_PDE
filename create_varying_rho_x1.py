# write some simple code to sketch pictures of my nondimensional model

# get imports
import numpy as np
import scipy.sparse as scp
# from matplotlib import pyplot as plt
from scipy.stats import norm, truncnorm, truncexpon, uniform
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Computer Modern"]})
from scipy.optimize import curve_fit

# function to make matrices for C
def make_matrix_C(delta_x, n_rows, n_cols,D):
    main_diag = (-4/delta_x**2)*np.ones((n_rows,))
    
    above_diag = np.empty((n_rows-1,))
    above_diag[0] = 2/delta_x**2 
    above_diag[1:] = 1/delta_x**2

    below_diag = np.empty((n_rows-1,))
    below_diag[-1] = 2/delta_x**2 
    below_diag[:-1] = 1/delta_x**2
    
    C_mat = D*scp.diags((below_diag, main_diag, above_diag), (-1,0,1))

    # kron products to obtain full matrix
    struct_main_diag = scp.eye(n_cols)
    main_diag_mat = scp.kron(struct_main_diag, C_mat)

    # do the off diagonal matrices 
    C_above_below = (D/delta_x**2)*scp.eye(n_rows)
    C_top_bottom = (2*D/delta_x**2)*scp.eye(n_rows)

    # make the structure matrices for these 
    above_diag_struct = np.empty(n_cols - 1)
    above_diag_struct[0] = 0
    above_diag_struct[1:] = 1
    below_diag_struct = np.empty(n_cols -1)
    below_diag_struct[-1] = 0
    below_diag_struct[:-1] = 1
    C_off_diag_struct = scp.diags((below_diag_struct, above_diag_struct),(-1,1))

    # top_bottom matrix
    top_struct = np.zeros(n_cols -1)
    top_struct[0] = 1
    bottom_struct = np.zeros(n_cols -1)
    bottom_struct[-1] = 1
    C_top_bottom_struct = scp.diags((top_struct, bottom_struct), (1,-1))

    # combine for full structure matrix for off diagonals
    mat_off_diags = scp.kron(C_off_diag_struct, C_above_below) + scp.kron(C_top_bottom_struct, C_top_bottom)

    # combine for full C matrix
    C_mat = main_diag_mat + mat_off_diags
    return C_mat

def make_matrix_C_advection(delta_x, n_rows, n_cols, D, adv):

    # make main diagonal matrix
    main_diag = D*(-4/delta_x**2)*np.ones((n_rows,))
    main_diag[1:-1] = main_diag[1:-1] - adv*(1/delta_x)
    
    above_diag = np.empty((n_rows-1,))
    above_diag[0] = D*2/delta_x**2 
    above_diag[1:] = D/delta_x**2

    below_diag = np.empty((n_rows-1,))
    below_diag[-1] = D*2/delta_x**2 
    below_diag[:-1] = D/delta_x**2
    below_diag[1:-1] = below_diag[1:-1] + adv*(1/delta_x)
    
    C_mat = scp.diags((below_diag, main_diag, above_diag), (-1,0,1))

    # kron products to obtain full matrix
    struct_main_diag = scp.eye(n_cols)
    main_diag_mat = scp.kron(struct_main_diag, C_mat)

    # do the off diagonal matrices 
    C_above_below = (D/delta_x**2)*scp.eye(n_rows)
    C_top_bottom = (2*D/delta_x**2)*scp.eye(n_rows)

    # make the structure matrices for these 
    above_diag_struct = np.empty(n_cols - 1)
    above_diag_struct[0] = 0
    above_diag_struct[0:] = 1
    below_diag_struct = np.empty(n_cols - 1)
    below_diag_struct[-1] = 0
    below_diag_struct[:-1] = 1
    C_off_diag_struct = scp.diags((below_diag_struct, above_diag_struct),(-1,1))

    # top_bottom matrix
    top_struct = np.zeros(n_cols -1)
    top_struct[0] = 1
    bottom_struct = np.zeros(n_cols -1)
    bottom_struct[-1] = 1
    C_top_bottom_struct = scp.diags((top_struct, bottom_struct), (-1,1))

    # combine for full structure matrix for off diagonals
    mat_off_diags = scp.kron(C_off_diag_struct, C_above_below) + scp.kron(C_top_bottom_struct, C_top_bottom)

    # combine for full C matrix
    C_mat = main_diag_mat + mat_off_diags
    return C_mat


# timestepping function for C
def timestep_C(A_vector, C_vector, delta_t, C_mat, gamma, beta):
    C_t_plus = C_vector + delta_t*(C_mat*C_vector - C_vector + gamma*A_vector**2/(beta**2 + A_vector**2))
    return C_t_plus

# timestepper for A
def timestep_A(A_vector, C_vector, delta_t, epsilon, rho_0, kappa):
    A_t_plus = A_vector + delta_t*(epsilon*C_vector + rho_0*A_vector - kappa*A_vector**2)
    return A_t_plus

# solve and plotter for output
def solve_pde(inits, n_time_points = 10000, domain_length = 10, delta_t = 1e-3, D = 1, gamma = 1, beta = 1, epsilon = 1e-2, rho_0 = 1.2, kappa = 1, adv = 0):

    # set delta_x
    n_rows = len(inits[0][:,0])
    n_cols = len(inits[0][0,:])
    delta_x = domain_length/(n_rows - 1)

    # create C matrix
    if adv != 0:
        C_mat = make_matrix_C_advection(delta_x, n_rows, n_cols, D, adv)
    else:
        C_mat = make_matrix_C(delta_x, n_rows, n_cols, D)
    
    # lambda-ize A and C functions
    lambda_timestepper = lambda x: (timestep_A(x[0],x[1], delta_t, epsilon, rho_0, kappa), timestep_C(x[0],x[1], delta_t, C_mat, gamma, beta))

    # start the loop
    cur_soln = [np.reshape(x, (n_rows*n_cols, 1)) for x in inits]
    soln_list = [cur_soln]
    for time_count in range(n_time_points):
        cur_soln = lambda_timestepper(cur_soln)
        soln_list += [cur_soln]
    

    plot_array_A = soln_list[-1][0]

    return np.reshape(plot_array_A,(n_rows*n_cols,1)), soln_list

# create a field of varying parameter values
def create_varying_parameter(inits, mean = 1, std_dev = 0.1):
    a, b = (0.1 - mean) / std_dev, (1.9 - mean) / std_dev
    non_smoothed_noise = truncnorm.rvs(a,b,loc = mean, scale = std_dev, size = (inits.shape[0] + 2, inits.shape[1] +2))
    out = np.empty_like(non_smoothed_noise)
    for i in range(1,inits.shape[0]+1):
        for j in range(1, inits.shape[1]+1):
            out[i,j] = np.mean(non_smoothed_noise[i-1:i+2,j-1:j+2])
    out = out[1:-1,1:-1]
    # plt.imshow(out)
    # plt.show()
    return np.reshape(out,(out.shape[0]*out.shape[1],1))

# get the mean and variance of the distance of invaded cells from the origin as a function of 
# time.
def get_var_invaded(soln_list, threshold = 0.5):
    x_vals = np.linspace(-5,5,101)
    y_vals = np.linspace(-5,5,101)
    xv, yv = np.meshgrid(x_vals, y_vals)
    xv_1D = np.reshape(xv, (101**2,1))
    yv_1D = np.reshape(yv, (101**2,1))
    var_list = []
    mean_list = []
    for soln in soln_list:
        cur_soln_1d = np.reshape(soln[0], (101**2,1))
        invaded_cells_x = xv_1D[cur_soln_1d > 0.5]
        invaded_cells_y = yv_1D[cur_soln_1d > 0.5]

        dist_cells = np.sqrt(invaded_cells_x**2)
        var_list += [np.var(dist_cells)]
        mean_list += [np.mean(dist_cells)]
    
    return var_list, mean_list

# piecewise function used for data
def fit_func(x, slope, intercept):
    return np.reshape(np.maximum(0, slope*x + intercept),(10001,))

# main loop
if __name__ == "__main__":
    # Create the initial condition
    inits1_A = np.zeros((101,101))
    inits1_A[:,50:51] = 1
    inits_1_C = np.zeros((101,101))

    # do 1000 iterates of the model
    max_its = 1000
    mean_list_out = np.empty((10001,max_its))
    threshold_mean = np.empty((9,max_its))

    for its in range(max_its):
        gamma_mod = create_varying_parameter(inits1_A)
        one_d_out, soln_list = solve_pde((inits1_A,inits_1_C), D = 1, adv = 0, rho_0 = gamma_mod, kappa = gamma_mod)

        x_vals = np.linspace(-5,5,101)
        y_vals = np.linspace(-5,5,101)
        xv, yv = np.meshgrid(x_vals, y_vals)
        xv_1D = np.reshape(xv, (101**2,1))
        yv_1D = np.reshape(yv, (101**2,1))

        invaded_cells_x = xv_1D[one_d_out > 0.5]
        invaded_cells_y = yv_1D[one_d_out > 0.5]

        dist_cells = np.sqrt(invaded_cells_x**2)

        var_list, mean_list = get_var_invaded(soln_list)
        mean_list_out[:,its] = mean_list

        invaded_per_threshold = []
        for decile in range(1,10):
            invaded_per_threshold += [len(one_d_out[one_d_out > decile*0.1])/len(one_d_out)]
        threshold_mean[:,its] = invaded_per_threshold

    # average across all 1000 entries
    mean_mean_list = np.mean(mean_list_out, axis = 1)
    mean_threshold = np.mean(threshold_mean, axis = 1)
    t_vals = np.reshape(np.linspace(0,10,10001), (10001, 1))
    curve, covar = curve_fit(fit_func, t_vals, mean_mean_list, p0=(0.35,-3))
    # print coefficients of the piecewise function
    print(curve)

    # fig1, (ax1,ax2,ax3) = plt.subplots(3,1)
    # ax1.imshow(np.reshape(one_d_out,(101,101)), origin = 'lower', aspect = 'auto')

    np.savetxt('mean_mean_list_rho_x1.csv', np.concatenate((np.reshape(np.linspace(0,10,10001), (10001, 1)), np.reshape(np.array(mean_mean_list), (10001,1))), axis = 1))
    np.savetxt('A10_rep_rho_x1.csv', np.concatenate((xv_1D, yv_1D, one_d_out), axis = 1))
    np.savetxt('invaded_per_threshold_rho_x1.csv', np.concatenate((np.reshape(np.array(mean_threshold), (9,1)), 0.1*np.reshape(np.array(range(1,10)), (9,1))), axis = 1))

    # ax2.plot(t_short, fit)
    # ax2.plot(t_short, mean_mean_list_short)

    # ax3.plot(0.1*np.array(range(1,10)), mean_threshold)
    # plt.show()

    
    




