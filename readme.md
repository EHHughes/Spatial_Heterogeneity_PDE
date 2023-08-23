This repository contains the data used in Hughes et al (under review). Representative code used to generate the data is also provided. Please contact Elliott Hughes at elliott.hughes@rebuen.ox.ac.uk for more information.

<font size = 5>Index<font size = 2>

**create_varying_rho_x1.py**

This file creates the data for $\rho_0 \sim S(X_1)$ seen in figures 12 and 13. With appropriate modifications (e.g. changing $\rho_0$ to $\kappa$, etc) it can be used to create the data for figures 1-17 except where $\Delta_t = 10^{-2}$. Note that in general results will vary slightly from those shown in the paper, as no random seed is provided.

**fine_t.py**

This file handles the case where $\Delta_t = 10^{-2}$ in figure 8. Otherwise it is essentially identical to the previous file.

**A10_rep_const.csv**

Corresponds to the representative state given in figure 1. The columns in these .csv's are values of x,y, and A respectively.

**Gamma_Output.csv**

Corresponds to figure 3 (this csv is formatted as above).

**A10_rep.csv**

Corresponds to figure 4 (see above for formatting).

**A10_rep_gamma/kappa/rho.csv**

These files correspond to the representative states discussed in figures 
13, 15, and 17 (see above for formatting).

**invaded_per_threshold\*.csv**

These files correspond to figures 5, 12, 14, and 16. The first column gives the proportion of cells invaded at $t=10$ for varying $T$ while the second column gives the corresponding values of $T$.

**dist_cells\*.csv**

Corresponds to figure 6. Each entry is the distance between an invaded cell and the initial condition.

**mean\*\_list_\*\.csv**

Corresponds to figure 7-12, 14, and 16. The 
first column gives the average distance of 
invaded cells to the initial condition (or the curve fitted to this data for mean_list_fitted_curve_const_rho.csv or 
mean_list_fitted_curve_vary_rho.csv), 
while the second column gives the current time.




