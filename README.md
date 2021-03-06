# 2D_Proliferation-Invasion-Recruitment
These MATLAB codes were used to generate the simulations in [Massey, S.C., Rockne, R.C., Hawkins-Daarud, A., Gallaher, J., Anderson, A.R.A, Canoll, P., Swanson, K.R. Simulating PDGF-driven Glioma Growth and Invasion in an Anatomically Accurate Brain Domain, *Bull Math Biol*, 80(5):1292-1309, 2018](https://doi.org/10.1007/s11538-017-0312-3).

The mathematical model is solved using Godunov splitting to apply the diffusion and reaction operators. We also split to set up the problem as locally one-dimensional, solving in x- and y-directional sweeps.

File run_2dpir_sims_set_newR0.m is the code where the simulation parameters are set, and is used to call the primary simulation function, run_pir_sim_newR0.m. This function then calls the other codes, set_interfaces.m and pir_numerics_2d.m

set_interfaces.m helps set up the locally one dimensional approach to solve in x-sweeps and y-sweeps.

pir_numerics_2d.m is the main numerical solving code and runs the diffusion solves, calling diffusion_op.m, as well as the reaction solves, calling pir_reaction.m 

At the end of each time step, pir_numerics_2d.m calls T2radius.m to compute the radius of the area of elevated cellular density that would be detected on T1 and T2 MRI modalities, to determine if this is a size at which we would like to save simulation solutions (if saving by size - if saving by time, this radius is recorded at the time points at which the simulation saves solution data).
