Notes:

I have implemented two different interior point methods:

	1) Primal-Dual Infeasible Start Interior Point Method (PD IIPM)

	2) Primal Barrier Infeasible Start Interior Point Method (PB IIPM).


If possible, the PD IIPM should be used because 
	a) much faster 
	b) no PhaseI() needed
	c) easier to handle (fewer parameters to be specified)
	
-George