# write some simple code to sketch pictures of my nondimensional model

# get imports
import numpy as np
import scipy.sparse as scp
from scipy.optimize import curve_fit
# from matplotlib import pyplot as plt
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Computer Modern"]})

# function to make matrices for C
def make_matrix_C(delta_x, n_pts,D):
    main_diag = (-2/delta_x**2)*np.ones((n_pts,))
    
    above_diag = np.empty((n_pts-1,))
    above_diag[0] = 2/delta_x**2 
    above_diag[1:] = 1/delta_x**2

    below_diag = np.empty((n_pts-1,))
    below_diag[-1] = 2/delta_x**2 
    below_diag[:-1] = 1/delta_x**2
    
    C_mat = D*scp.diags((below_diag, main_diag, above_diag), (-1,0,1))
    return C_mat

def make_matrix_C_advection(delta_x, n_pts, D, adv):
    main_diag = D*(-2/delta_x**2)*np.ones((n_pts,))
    main_diag[1:-1] = main_diag[1:-1] - adv*(1/delta_x)
    
    above_diag = np.empty((n_pts-1,))
    above_diag[0] = D*2/delta_x**2 
    above_diag[1:] = D/delta_x**2

    below_diag = np.empty((n_pts-1,))
    below_diag[-1] = D*2/delta_x**2 
    below_diag[:-1] = D/delta_x**2
    below_diag[1:-1] = below_diag[1:-1] + adv*(1/delta_x)
    
    C_mat = scp.diags((below_diag, main_diag, above_diag), (-1,0,1))
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
def solve_pde(inits, n_time_points = 999999, domain_length = 10, delta_t = 1e-5, D = 1, gamma = 1, beta = 1, epsilon = 1e-2, rho_0 = 1.2, kappa = 1, adv = 0):

    # set delta_x
    n_pts = len(inits[0])
    delta_x = domain_length/(n_pts - 1)

    # create C matrix
    if adv != 0:
        C_mat = make_matrix_C_advection(delta_x, n_pts, D, adv)
    else:
        C_mat = make_matrix_C(delta_x, n_pts, D)
    
    # lambda-ize A and C functions
    lambda_timestepper = lambda x: (timestep_A(x[0],x[1], delta_t, epsilon, rho_0, kappa), timestep_C(x[0],x[1], delta_t, C_mat, gamma, beta))

    # start the loop
    cur_soln = inits
    soln_list = [inits]
    for time_count in range(n_time_points):
        cur_soln = lambda_timestepper(cur_soln)
        soln_list += [cur_soln]
    
    # plot the solutions
    plot_array = np.array(soln_list)
    plot_array_A = plot_array[:,0,:]
    plot_array_C = plot_array[:,1,:]

    return plot_array_A

# get the evolution of the mean and variance of distances to the initial 
# condition over time.
def get_var_invaded(numerical_soln, x_vals, threshold = 0.5):
    var_list = []
    mean_list = []
    for i in range(numerical_soln.shape[0]):
        soln = numerical_soln[i,:]
        invaded_cells_x = x_vals[soln > 0.5]

        dist_cells = np.sqrt(invaded_cells_x**2) - 0.05
        dist_cells[dist_cells < 0] = 0
        var_list += [np.var(dist_cells)]
        mean_list += [np.mean(dist_cells)]
    
    return var_list, mean_list

def fit_func(x, slope, intercept):
    return np.reshape(np.maximum(0, slope*x + intercept),(1000000,))

# main loop
if __name__ == "__main__":

    inits1_A = np.zeros((1001,))
    inits1_A[495:506] = 1
    inits_1_C = np.zeros((1001,))

    numerical_soln = solve_pde((inits1_A,inits_1_C), D = 1, adv = 0, rho_0 = 1)

    x_vals = np.linspace(-5,5,1001)

    invaded_cells_x = x_vals[numerical_soln[-1,:] > 0.5]

    dist_cells = np.sqrt(invaded_cells_x**2)

    var_list, mean_list = get_var_invaded(numerical_soln, x_vals)

    invaded_per_threshold = []
    for decile in range(1,10):
        invaded_per_threshold += [len(numerical_soln[-1,:][numerical_soln[-1,:] > decile*0.1])/len(numerical_soln[-1,:])]


    t_vals = np.reshape(np.linspace(0,10,1000000), (1000000, 1))
    mean_array = np.array(mean_list)
    curve, covar = curve_fit(fit_func, t_vals, mean_list, p0=(0.35,-3))
    print(curve)

    np.savetxt('mean_list_t_fine.csv', np.concatenate((np.reshape(np.linspace(0,10,1000000), (1000000, 1)), np.reshape(mean_array, (1000000,1)), np.reshape(fit_func(t_vals, curve[0], curve[1]), (1000000,1))), axis = 1))

    var_array = np.array(var_list)

    np.savetxt('var_list_t_fine.csv', np.concatenate((np.reshape(np.linspace(0,10,1000000), (1000000, 1)), np.reshape(var_array, (1000000,1))), axis = 1))

    


