# Summary

This repository contains all files needed to reproduce the paper "Computing Absorbing Times via Fluid Approximations", accepted Advances in Applied Probability (to appear in september 2017). 

# Recipe

To generate the figures and the paper, type: 
```{sh}
make
```

*Warning* : to fasten the above process, some simulation results are stored in the files simu/dict_dist_N_K_coupon and simu/dict_dist_N_K_in_a_row

If you want to erase these files and regenerate all simulations, type :
```{sh}
make redo_all_simu
```
(this should take about 2h)
