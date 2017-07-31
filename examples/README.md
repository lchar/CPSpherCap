These are the Closest point simulation files that we use

List of files
------------------------------------------------
      SpherCap_CrossSection.m
	2D Closest point on a full cross-section (Solutions are currently not symmetric)
	
      SpherCap_EvolvingDomain.m
	2D Closest point on a half cross-section (there is currently a bug)
	
      SpherCap_EvolvingDomainXY.m
	2D Closest point method with XY system for spherical cap 
	[This will probably get abandonned in favour of the 3D approach of the next three files]
	
      example_RD_cap_Bruss_Dir.m
	3D Closest point method with explicit time-stepping
	
      example_RD_cap_Bruss_Dir_IMEX.m
	3D Closest point method with implicit time-stepping (for the laplacian)
	This one currently works and will be adapted to the Evolving domain problem.
	
      RD_Cap_EvolvingDomain.m 	             
	3D Closest point method with explicit time-stepping and evolving domain
