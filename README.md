# wedge
Fortran codes related to the study of two dimensional quantum scattering inside a wedge waveguide via Boundary Wall Method

gaintwocircs.f90  
Objective: Write in a file the values of the circles parameters, and the respecive mean values related to the segments of the circles.

wedgearbresv.f90
Objective: Write in a file the mean value of the squared module of the T-matrix elements for each k. This values are used to find the numeric resonant k close to the initial guess "kanalitico". We vary the search range of k to meet the desired precision.

wedgeincv.f90
Objective: Write in a file the wavefunction scattered by a hard wall
