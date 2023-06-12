pro Example_stdZeta_noDEM_SingleThread
 ;Test code: computing the gyroresonance and free-free emission using pre-defined zeta function and isothermal plasma, 
 ;for a single line-of-sight
 
 libname='GRFF_DEM_Transfer_64.dll' ;library name (change if necessary)
 
 Nz=100L ;number of voxels
 Nf=240L ;number of frequencies
 
 ;array of dimensions etc.:
 Lparms=[Nz,    $ ;number of voxels
         Nf,    $ ;number of frequencies
         0,     $ ;number of temperatures (0 - no DEM/DDM)
         1,     $ ;global DEM on/off key (off)
         1      ] ;global DDM on/off key (off)
 
 ;global floating-point parameters: 
 Rparms=[1d18,  $ ;source area, cm^2
         3d9,   $ ;starting frequency, Hz 
         0.005  ] ;logarithmic step in frequency 

 ;single-voxel parameters: 
 ParmLocal=dblarr(15)
 ParmLocal[0]=2d9    ;source depth, cm (total depth - the depths for individual voxels will be computed later)
 ParmLocal[1]=3d6    ;plasma temperature, K
 ParmLocal[2]=9d8    ;electron concentration, cm^{-3}
 ParmLocal[3]=100d0  ;magnetic field, G (will be changed later)
 ParmLocal[4]=120d0  ;viewing angle, degrees
 ParmLocal[5]=0d0    ;azimuthal angle, degrees
 ParmLocal[6]=0      ;emission mechanism flag (all on)
 ParmLocal[7]=30     ;maximum harmonic number
 ParmLocal[8]=0      ;proton concentration, cm^{-3} (not used in this example)
 ParmLocal[9]=0      ;neutral hydrogen concentration, cm^{-3}
 ParmLocal[10]=0     ;neutral helium concentration, cm^{-3}
 ParmLocal[11]=1     ;local DEM on/off key (off)
 ParmLocal[12]=1     ;local DDM on/off key (off)
 ParmLocal[13]=0     ;element abundance code (coronal, following Feldman 1992)
 ParmLocal[14]=0     ;reserved
 
 ;multi-voxel parameters along line-of-sight:
 Parms=dblarr(15, Nz)
 for i=0, Nz-1 do begin
  Parms[*, i]=ParmLocal
  
  Parms[0, i]=ParmLocal[0]/Nz ;depths of individual voxels
  Parms[3, i]=1000d0-700d0*i/(Nz-1) ;magnetic field decreases linearly from 1000 to 300 G
  
  ;all other parameters are the same in all voxels
 endfor
 
 RL=dblarr(7, Nf) ;output array
 
 ;main call:
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Parms, 0, 0, 0, RL)
 ;temperature, DEM and DDM arrays are not used, so the corresponding parameters can be set to any value

 ;extracting the results:
 freq=reform(RL[0, *]) ;frequencies
 I_L=reform(RL[5, *]) ;left-hand-polarized flux, exact mode coupling
 I_R=reform(RL[6, *]) ;right-hand-polarized flux, exact mode coupling
 
 ;--------------------------------------------
 ;Now plotting the results
 
 window, 1, title='Intensity'
 wset, 1
 plot, freq, I_L+I_R, /xlog, /ylog, xstyle=1, $
       xtitle='Frequency, GHz', ytitle='Intensity, sfu'
 
 window, 2, title='Polarization degree'
 wset, 2
 plot, freq, (I_L-I_R)/(I_L+I_R), /xlog, yrange=[-1, 1], xstyle=1, ystyle=1, $
       xtitle='Frequency, GHz', ytitle='Polarization degree'
end
