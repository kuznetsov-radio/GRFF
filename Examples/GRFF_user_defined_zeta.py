import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav

def initGET_MW(libname):
	"""
	Python wrapper for fast gyrosynchrotron codes.
	Identical to GScodes.py in https://github.com/kuznetsov-radio/gyrosynchrotron
	This is for the single thread version
	@param libname: path for locating compiled shared library
	@return: An executable for calling the GS codes in single thread
	"""
	_intp = ndpointer(dtype=ctypes.c_int32, flags='F')
	_doublep = ndpointer(dtype=ctypes.c_double, flags='F')

	libc_mw = ctypes.CDLL(libname)
	mwfunc = libc_mw.PyGET_MW_USER  ### use PyGET_MW_USER for user defined zeta values
	mwfunc.argtypes = [_intp, _doublep, _doublep, _doublep, _doublep,_doublep,_doublep,_doublep,_doublep,_doublep]
	mwfunc.restype = ctypes.c_int
	return mwfunc
    


libname='/home/surajit/GR-FF_general/Binaries/GRFF_DEM_Transfer.so'
abundance_file='/home/surajit/GR-FF_general/Examples/zeta_Feldman.sav'
Nf=100    # number of frequencies
NSteps=1  # number of nodes along the line-of-sight
N_temp_DEM=2
Npix=1
T=np.array([5*1e5,1e6],dtype='double')


GET_MW = initGET_MW(libname)
abundance_data=readsav(abundance_file)

fzeta_arr=np.asfortranarray(abundance_data['frqhz'],dtype='double')
Tzeta_arr=np.asfortranarray(abundance_data['T'],dtype='double')
zeta_arr=np.asfortranarray(abundance_data['eta'],dtype='double')

Nf_zeta=np.size(fzeta_arr) #number of frequencies
NT_zeta=np.size(Tzeta_arr) #number of temperatures
N_zeta=1
 
Lparms=np.zeros(8, dtype='int32') # array of dimensions etc.
Lparms[0]=NSteps
Lparms[1]=Nf
Lparms[2]=N_temp_DEM
Lparms[3]=1
Lparms[4]=1
Lparms[5]=Nf_zeta
Lparms[6]=NT_zeta
Lparms[7]=N_zeta

 
Rparms=np.zeros(5, dtype='double') # array of global floating-point parameters - for a single LOS
Rparms[0]=1e20 # area, cm^2
Rparms[1]=1e9  # starting frequency to calculate spectrum, Hz
Rparms[2]=0.02 # logarithmic step in frequency
Rparms[3]=0   # f^C
Rparms[4]=0   # f^WH

depth_cm=1e10

ParmLocal = np.zeros(24, dtype='double')  # array of voxel parameters - for a single voxel
ParmLocal[0] = depth_cm / NSteps  # voxel depth, cm
ParmLocal[1] = 2*1e6  # T_0, K
ParmLocal[2] = 1e10  # n_0 - thermal electron density, cm^{-3}
ParmLocal[3] = 200  # B - magnetic field, G
ParmLocal[6] = 0  # distribution over energy (PLW is chosen)
ParmLocal[7] = 1e5  # n_b - nonthermal electron density, cm^{-3}
ParmLocal[9] = 10*1e-3  # E_min, MeV
ParmLocal[10] = 10  # E_max, MeV
ParmLocal[12] = 4  # \delta_1
ParmLocal[14] = 0  # distribution over pitch-angle (isotropic is chosen)
ParmLocal[15] = 90  # loss-cone boundary, degrees
ParmLocal[16] = 0.2  # \Delta\mu

Parms = np.zeros((24, NSteps), dtype='double', order='F')  # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i] = ParmLocal  # most of the parameters are the same in all voxels
    
RL = np.zeros((7, Nf), dtype='double', order='F')  # input/output array

DEM_arr=np.zeros((N_temp_DEM,NSteps),dtype='double',order='F')
DDM_arr=np.zeros((N_temp_DEM,NSteps),dtype='double',order='F')
DDM_arr[...]=DEM_arr[...]


res = GET_MW(Lparms, Rparms, Parms,T,DEM_arr,DDM_arr,fzeta_arr, Tzeta_arr,zeta_arr, RL)

RR=RL[5]
LL=RL[6]
freqs=RL[0]

I=RR+LL

plt.loglog(freqs,I)
plt.show()

  
