# write some simple code to sketch pictures of my nondimensional model

# get imports
import numpy as np
import scipy.sparse as scp
from matplotlib import pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern"]})

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
def solve_pde(inits, n_time_points = 19999, domain_length = 20, delta_t = 1e-3, D = 1, gamma = 1, beta = 1, epsilon = 1e-2, rho_0 = 1.2, kappa = 1, adv = 0):

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
    fig, (ax1, ax2) = plt.subplots(2,1)
    plot_array = np.array(soln_list)
    plot_array_A = plot_array[:,0,:]
    plot_array_C = plot_array[:,1,:]
    kym1 = ax1.imshow(plot_array_A, aspect = 'auto', origin = 'lower')
    #ax1.set_xticks([0,25,50,75,100],[-5,-2.5,0,2.5,5])
    # ax1.set_yticks([0,5000,10000,15000,20000,25000,30000],[0,50,100,150,200,250,300])
    ax1.set_xlabel('$x$', loc = 'right', fontsize = 18)
    ax1.set_ylabel('Years', loc = 'top', rotation = 'horizontal', fontsize = 18)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    plt.colorbar(kym1, ax = ax1)

    kym2 = ax2.imshow(plot_array_C, aspect = 'auto', origin = 'lower')
    #ax2.set_xticks([0,25,50,75,100],[-5,-2.5,0,2.5,5])
    # ax2.set_yticks([0,5000,10000,15000,20000,25000,30000],[0,50,100,150,200,250,300])
    ax2.set_xlabel('$x$', loc = 'right', fontsize = 18)
    ax2.set_ylabel('Years', loc = 'top', rotation = 'horizontal', fontsize = 18)
    ax2.tick_params(axis = 'x', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    plt.colorbar(kym2, ax = ax2)

    plt.tight_layout()
    

    plt.show()

    return plot_array_A


# main loop
if __name__ == "__main__":
    rho_choice = 1
    inits1_A = np.zeros((201,))
    inits1_A[100:101] = rho_choice
    inits_1_C = np.zeros((201,))

    x_vals = np.linspace(-10,10,201)
    y_vals = np.linspace(0,20,200)
    xv, yv = np.meshgrid(x_vals, y_vals)

    numerical_soln = solve_pde((inits1_A,inits_1_C), D = 1, adv = 0, rho_0 = rho_choice)
    numerical_soln = numerical_soln[::100,:]
    out = np.concatenate((np.reshape(xv, (numerical_soln.shape[0]*numerical_soln.shape[1],1)),np.reshape(yv, (numerical_soln.shape[0]*numerical_soln.shape[1],1)),np.reshape(numerical_soln, (numerical_soln.shape[0]*numerical_soln.shape[1],1))), axis = 1)
    np.savetxt("num_soln_travelling_wv.csv",out)

    