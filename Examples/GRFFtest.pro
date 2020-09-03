pro GRFFtest
 Nz=100L ;number of voxels

 restore, 'DEMexample.sav' ;reading the DEM, DDM and T arrays - for a single voxel
 T_arr=double(T_arr)
 NT=n_elements(T_arr)
 ;preparing the NT*Nz arrays (with the same DEM and DDM in all voxels):
 DEM_arr=double(rebin(DEM_arr, NT, Nz))
 DDM_arr=double(rebin(DDM_arr, NT, Nz))

 ;single-voxel parameters: 
 Parms=dblarr(15)
 Parms[0]=1d18   ;source area, cm^2
 Parms[1]=2d9    ;source depth, cm (total depth - the depths for individual voxels will be computed later)
 Parms[2]=1d6    ;plasma temperature, K (not used in this example, because DDM-based temperature is used)
 Parms[3]=1d9    ;electron density, cm^{-3} (not used in this example, because DDM-based density is used)
 Parms[4]=100d0  ;magnetic field, G (will be changed later)
 Parms[5]=120d0  ;viewing angle, degrees
 Parms[6]=0d0    ;azimuthal angle, degrees
 Parms[7]=3d9    ;starting frequency, Hz 
 Parms[8]=0.005  ;logarithmic step in frequency 
 Parms[9]=0      ;emission mechanism flag (all on)
 Parms[10]=30    ;maximum harmonic number
 Parms[11]=0     ;neutral hydrogen concentration, cm^{-3}
 Parms[12]=0     ;neutral helium concentration, cm^{-3}
 Parms[13]=0     ;DEM on/off key (on)
 Parms[14]=0     ;abundance (coronal)
 
 ;multi-voxel parameters along line-of-sight:
 ParmsM=dblarr(15, Nz)
 for i=0, Nz-1 do begin
  ParmsM[*, i]=Parms
  
  ParmsM[1, i]=Parms[1]/Nz ;depths of individual voxels
  ParmsM[4, i]=1000d0-700d0*i/(Nz-1) ;magnetic field decreases linearly from 1000 to 300 G
  
  ;all other parameters are the same in all voxels
 endfor
 
 Nf=240L ;number of frequencies
 
 Ndat=[Nz, Nf, NT, n_elements(Parms)] ;array of dimensions
 
 RL=dblarr(7, Nf) ;output array
 
 libname='MW_Transfer64.dll'
 
 ;main call:
 res=call_external(libname, 'GET_MW', Ndat, ParmsM, T_arr, DEM_arr, DDM_arr, RL, /unload)

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