This repository contains the data used in Hughes et al (under review). Representative code used to generate the data is also provided. Please contact Elliott Hughes at elliott.hughes@maths.ox.ac.uk for more information.

<font size = 5>Index<font size = 2>

**create_varying_rho_x1.py**

This file creates the data for $\rho_0 \sim S(X_1)$ seen in figures 11 and 12. With appropriate modifications (e.g. changing $\rho_0$ to $\kappa$, etc) it can be used to create the data for figures 1-16 except where $\Delta_t = 10^{-2}$. Note that in general results will vary slightly from those shown in the paper, as no random seed is provided.

**fine_t.py**

This file handles the case where $\Delta_t = 10^{-2}$ in figure 8. Otherwise it is essentially identical to the previous file.

**advective_teardrop_vary_rho_0.py**

Produces the output for figure 17.

**toy_sim_double_domain**

Code to produce the 1D solution on a double width domain (e.g. figure S1). With appropriate modifications it can also be used to reproduce figure 1.

**Numerical_Output.csv**

Data for the output in figure 1. The coluns in this .csv are x,t and A respectively.

**A10_rep_const.csv**

Corresponds to the representative state given in figure 2. The columns in these .csv's are values of x,y, and A respectively.

**Random_Field_Test_Output.csv**

Corresponds to figure 4 (this csv is formatted as above).

**A10_rep_fig_5.csv**

Corresponds to figure 5 (right hand figure, see above for formatting).

**A7_rep_fig_5.csv**

Corresponds to figure 5 (left hand figure, see above for formatting).

**A10_rep_gamma/kappa/rho.csv**

These files correspond to the representative states discussed in figures 12, 14, and 16 (see above for formatting).

**invaded_per_threshold\*.csv**

These files correspond to figures 6, 11, 13, and 15. The first column gives the proportion of cells invaded at $t=10$ for varying $T$ while the second column gives the corresponding values of $T$.

**dist_cells.csv** and **dist_cells_flat.csv**

Corresponds to figure 7. Each entry is the distance between an invaded cell and the initial condition for varying rho_0 (dist_cells.csv) and for constant rho_0 (dist_cells_flat.csv).

**mean\*\_list_\*\.csv**

Corresponds to figure 8-10, 11, 13, and 15. The first column gives the average distance of invaded cells to the initial condition (or the curve fitted to this data for mean_list_fitted_curve_const_rho.csv or mean_list_fitted_curve_vary_rho.csv), while the second column gives the current time.

**mean_list_t_fine.csv**

Fine t (e.g. $\Delta_t = 10^{-2}$) output for figures 8 and 9 (see fine_t.py for code).

**upper/lower\*_int_list\*.csv**

Gives the upper and lower intervals for figure 10 (upper/lower_init_list.csv) and S2-S10. Uses the same format as mean_list.csv.

**advection_with_varying_rho0.csv**

Gives the output for figure 17. Uses the same format as A10_rep_gamma_x1.csv.

**num_soln_travelling_wv.csv**

Gives the evolution of the 1D state for figure S1. Uses the same format as figure 1.




