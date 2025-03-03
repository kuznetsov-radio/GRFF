#pragma once

int InterpolateEBTEL(int NQ, int NL, int NT, double Q, double L, 
	                 float *Qrun, float *Lrun, float *DEM_run, float *DDM_run, 
	                 double *DEM, double *DDM);