"inlist_extra" - change mass and/or metallicity here
"inlist_project" - main inlist
	- change x_ctrl(18) under controls for efficiency factor for eruptive mass loss
	- Ensure use_other_wind = .true. when using eruptive mass loss model
	- Ensure both use_superad_reduction = .false. and okay_to_reduce_gradT_excess = .false. when using eruptive mass loss model
"src/run_star_extras_dejager.f90" - run_stars_extras without decin low-T mass loss change