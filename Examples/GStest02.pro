pro GStest02
 ;Test code: homogeneous source, free-free, neutrals, Saha equation, no DEM, no DDM. no GR

 Nz=1L ;number of voxels
 Nf=150L ;number of frequencies
 
 ;array of dimensions etc.:
 Lparms=[Nz,    $ ;number of voxels
         Nf,    $ ;number of frequencies
         0,     $ ;number of temperatures
         1,     $ ;DEM on/off key (off)
         1      ] ;DDM on/off key (off)
 
 ;global floating-point parameters: 
 Rparms=[1d20,  $ ;source area, cm^2
         1d9,   $ ;starting frequency, Hz 
         0.02  ] ;logarithmic step in frequency 

 ;single-voxel parameters: 
 Parms=dblarr(15)
 Parms[0]=1d7    ;source depth, cm
 Parms[1]=5d3    ;plasma temperature, K
 Parms[2]=1d13   ;total atomic electron density, cm^{-3}
 Parms[3]=0.2    ;magnetic field, G (will be changed later)
 Parms[4]=60d0   ;viewing angle, degrees
 Parms[5]=0d0    ;azimuthal angle, degrees
 Parms[6]=1      ;emission mechanism flag
 Parms[7]=0      ;maximum harmonic number
 Parms[8]=0      ;proton concentration, cm^{-3}
 Parms[9]=0      ;neutral hydrogen concentration, cm^{-3}
 Parms[10]=0     ;neutral helium concentration, cm^{-3}
 Parms[11]=1     ;DEM on/off key (off)
 Parms[12]=1     ;DDM on/off key (off)
 Parms[13]=1     ;abundance (chromospheric Caffau)
 Parms[14]=0     ;reserved
 
 RL=dblarr(7, Nf) ;output array
 
 libname='GRFF_DEM_Transfer_64.dll'
 
 RL=dblarr(7, Nf) ;output array
  
 ;main call:
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Parms, 0, 0, 0, RL, /unload)

 ;extracting the results:
 freq=reform(RL[0, *])
 I_L=reform(RL[5, *])
 I_R=reform(RL[6, *])

 ;--------------------------------------------------------------------
 
 window, 3, title='Intensity'
 wset, 3
 plot, freq, I_L+I_R, /xlog, /ylog, xstyle=1, yrange=[1d-2, 2d1], $
       xtitle='Frequency, GHz', ytitle='Intensity, sfu'
end