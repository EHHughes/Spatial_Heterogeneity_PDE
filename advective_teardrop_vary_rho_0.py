# write some simple code to sketch pictures of my nondimensional model

# get imports
import numpy as np
import scipy.sparse as scp
from matplotlib import pyplot as plt
from scipy.stats import truncnorm
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern"]})

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
    main_diag = (-4/delta_x**2)*np.ones((n_rows,)) 
    main_diag[1:-1] += -adv/delta_x
    
    above_diag = np.empty((n_rows-1,))
    above_diag[0] = 2/delta_x**2 
    above_diag[1:] = 1/delta_x**2

    below_diag = np.empty((n_rows-1,)) 
    below_diag[-1] = 2/delta_x**2 
    below_diag[:-1] = 1/delta_x**2 + adv/delta_x
    
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

def create_varying_parameter(inits, mean = 1, std_dev = 0.1):
    a, b = (0.1 - mean) / std_dev, (1.9 - mean) / std_dev
    non_smoothed_noise = truncnorm.rvs(a,b,loc = mean, scale = std_dev, size = (inits.shape[0] + 2, inits.shape[1] +2))
    out = np.empty_like(non_smoothed_noise)
    for i in range(1,inits.shape[0]+1):
        for j in range(1, inits.shape[1]+1):
            out[i,j] = 0.1+np.mean(non_smoothed_noise[i-1:i+2,j-1:j+2])
    out = out[1:-1,1:-1]
    # plt.imshow(out)
    # plt.show()
    return np.reshape(out,(out.shape[0]*out.shape[1],1))

# main loop
if __name__ == "__main__":
    x_vals = np.linspace(-5,5,101)
    y_vals = np.linspace(-5,5,101)
    xv, yv = np.meshgrid(x_vals, y_vals)
    xv_1D = np.reshape(xv, (101**2,1))
    yv_1D = np.reshape(yv, (101**2,1))

    inits1_blob_A = np.zeros((101,101))
    inits1_blob_A[48:52,47:53] = 1
    inits1_blob_A[47:53,48:52] = 1

    inits1_A = np.zeros((101,101))
    inits1_A[:,50:51] = 1
    inits_1_C = np.zeros((101,101))

    gamma_mod = create_varying_parameter(inits1_A, std_dev = 1)
    one_d_out, soln_list = solve_pde((inits1_blob_A,inits_1_C), D = 1, adv = 1, rho_0 = gamma_mod)

    fig1, ax1 = plt.subplots(1,1)
    ax1.imshow(np.reshape(one_d_out,(101,101)), origin = 'lower', aspect = 'auto')
    plt.show()
    np.savetxt('advection_with_varying_rho0.csv', np.concatenate((xv_1D, yv_1D, one_d_out), axis = 1))



    




