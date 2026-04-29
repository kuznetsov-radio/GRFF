import numpy as np
import matplotlib.pyplot as plt
import GRFFcodes # initialization library - located either in the current directory or in the system path

libname='./GRFF_DEM_Transfer_64.dll' # name of the executable library - located where Python can find it

GET_MW=GRFFcodes.initGET_MW(libname) # load the library

Nf=240     # number of frequencies
Nz=100     # number of nodes along the line-of-sight
 
Lparms=np.zeros(5, dtype='int32') # array of dimensions etc.
Lparms[0]=Nz
Lparms[1]=Nf
# other parameters are zero by default
 
Rparms=np.zeros(3, dtype='double') # array of global floating-point parameters
Rparms[0]=1e18  # area, cm^2
Rparms[1]=3e9   # starting frequency to calculate spectrum, Hz
Rparms[2]=0.005 # logarithmic step in frequency
 
ParmLocal=np.zeros(17, dtype='double') # array of voxel parameters - for a single voxel
ParmLocal[0]=2e7    # voxel size along the line-of-sight
ParmLocal[1]=3e6    #plasma temperature, K
ParmLocal[2]=9e8    #plasma density, cm^{-3}
#ParmLocal[3]=100.0  #magnetic field, G - will be set later
ParmLocal[4]=120.0  #viewing angle, degrees
ParmLocal[5]=0.0    #azimuthal angle, degrees
ParmLocal[6]=0      #emission mechanism flag (all on)
ParmLocal[7]=30     #maximum harmonic number
# other parameters are zero by default
 
Parms=np.zeros((17, Nz), dtype='double', order='F') # 2D array of input parameters - for multiple voxels
for i in range(Nz):
    Parms[:, i]=ParmLocal # most of the parameters are the same in all voxels
    Parms[3, i]=1000.0-700.0*i/(Nz-1) # magnetic field decreases linearly from 1000 to 300 G
 
RL=np.zeros((7, Nf), dtype='double', order='F') # input/output array
dummy=np.array(0, dtype='double')

# calculating the emission for analytical distribution (array -> off),
# the unused parameters can be set to any value
res=GET_MW(Lparms, Rparms, Parms, dummy, dummy, dummy, RL)
 
# retrieving the results
f=RL[0]
I_L=RL[5]
I_R=RL[6]

# plotting the results
plt.figure(1)
plt.plot(f, I_L+I_R)
plt.xscale('log')
plt.yscale('log')
plt.title('Total intensity')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')

r=1e-19*299792.458e5**2/(2.0*1.3806488e-16*(f*1e9)**2)/Rparms[0]*1.49599e13**2
Tb=(I_L+I_R)*r

plt.figure(2)
plt.plot(f, Tb)
plt.xscale('log')
plt.yscale('log')
plt.title('Brightness temperature')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')

plt.figure(3)
plt.plot(f, (I_L-I_R)/(I_L+I_R))
plt.xscale('log')
plt.title('Circular polarization degree')
plt.xlabel('Frequency, GHz')
plt.ylabel('Polarization degree')

plt.show()
