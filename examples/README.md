These are the Closest point simulation files that we use

List of files
------------------------------------------------
      eigf03.dat, eigf32.dat, eigf51.dat
	files containing patterns in band
	
      example_Bruss_sphercap_tor.m
	Brusselator on the spherical cap domain using the latest method
	
      example_RD_cap_Bruss_Dir_Evolving_SBDF_QP3.m
	Brusselator on the spherical cap domain using the old method (wide, single band skipping growth steps)
	
      example_RD_cap_Bruss_Dir.m
	3D Closest point method with explicit time-stepping
	
      example_RD_cap_Bruss_Dir_IMEX.m
	3D Closest point method with implicit time-stepping (for the laplacian)
	This one currently works and will be adapted to the Evolving domain problem.
	
      RD_Cap_EvolvingDomain.m 	             
	3D Closest point method with explicit time-stepping and evolving domain
