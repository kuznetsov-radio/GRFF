pro GRFFtest01
 ;Test code: gyroresonance, free-free, DEM, DDM, no neutrals, no Saha equation

 Nz=100L ;number of voxels
 Nf=240L ;number of frequencies
 
 restore, 'DEMexample.sav' ;reading the DEM, DDM and T arrays - for a single voxel
 T_arr=double(T_arr)
 NT=n_elements(T_arr)
 ;preparing the NT*Nz arrays (with the same DEM and DDM in all voxels):
 DEM_arr=double(rebin(DEM_arr, NT, Nz))
 DDM_arr=double(rebin(DDM_arr, NT, Nz))
 
 ;array of dimensions etc.:
 Lparms=[Nz,    $ ;number of voxels
         Nf,    $ ;number of frequencies
         NT,    $ ;number of temperatures
         0,     $ ;DEM on/off key (on)
         0      ] ;DDM on/off key (on)
 
 ;global floating-point parameters: 
 Rparms=[1d18,  $ ;source area, cm^2
         3d9,   $ ;starting frequency, Hz 
         0.005  ] ;logarithmic step in frequency 

 ;single-voxel parameters: 
 Parms=dblarr(15)
 Parms[0]=2d9    ;source depth, cm (total depth - the depths for individual voxels will be computed later)
 Parms[1]=1d6    ;plasma temperature, K
 Parms[2]=1d9    ;electron density, cm^{-3}
 Parms[3]=100d0  ;magnetic field, G (will be changed later)
 Parms[4]=120d0  ;viewing angle, degrees
 Parms[5]=0d0    ;azimuthal angle, degrees
 Parms[6]=0      ;emission mechanism flag
 Parms[7]=30     ;maximum harmonic number
 Parms[8]=0      ;total hydrogen concentration, cm^{-3}
 Parms[9]=0      ;neutral hydrogen concentration, cm^{-3}
 Parms[10]=0     ;neutral helium concentration, cm^{-3}
 Parms[11]=0     ;DEM on/off key (on)
 Parms[12]=0     ;DDM on/off key (on)
 Parms[13]=0     ;abundance (coronal)
 Parms[14]=0     ;reserved
 
 ;multi-voxel parameters along line-of-sight:
 ParmsM=dblarr(15, Nz)
 for i=0, Nz-1 do begin
  ParmsM[*, i]=Parms
  
  ParmsM[0, i]=Parms[0]/Nz ;depths of individual voxels
  ParmsM[3, i]=1000d0-700d0*i/(Nz-1) ;magnetic field decreases linearly from 1000 to 300 G
  
  ;all other parameters are the same in all voxels
 endfor
 
 RL=dblarr(7, Nf) ;output array
 
 libname='GRFF_DEM_Transfer_64.dll'
 
 ;main call:
 res=call_external(libname, 'GET_MW', Lparms, Rparms, ParmsM, T_arr, DEM_arr, DDM_arr, RL, /unload)

 ;extracting the results:
 freq=reform(RL[0, *])
 I_L=reform(RL[5, *])
 I_R=reform(RL[6, *])
 
 window, 1, title='Intensity'
 wset, 1
 plot, freq, I_L+I_R, /xlog, /ylog, xstyle=1, $
       xtitle='Frequency, GHz', ytitle='Intensity, sfu'
 
 window, 2, title='Polarization degree'
 wset, 2
 plot, freq, (I_L-I_R)/(I_L+I_R), /xlog, yrange=[-0.1, 1], xstyle=1, ystyle=1, $
       xtitle='Frequency, GHz', ytitle='Polarization degree'
end