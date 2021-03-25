The algorithms are coded in R

`corr_sim_slow.R` runs the Portkey algorithms along with the original Barker's (hence the name "slow").

`corr_shadow_sim.R` runs all algorithms including the shadow prior method for one run to produce the ACF plots.

`corr_optim.R` runs the Portkey and regular algorithms with approximately optimal scaling. This code will automatically quite in 24 hours. For the output, please look at the `corr_optim.Rout` file. 

`corr_plots.R` contains the codes that produce the plots.