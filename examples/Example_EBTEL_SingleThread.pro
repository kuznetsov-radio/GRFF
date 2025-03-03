pro Example_EBTEL_SingleThread
 libname='GRFF_DEM_Transfer_64.dll' 
 
 Nz=100L ;number of voxels
 Nf=240L ;number of frequencies
 
 ;load the EBTEL table:
 restore, gx_findfile('ebtel_ss.sav')
 T_arr=10d0^logtdem
 s=size(DEM_cor_run, /dimensions)
 NT=s[0]
 NQ=s[1]
 NL=s[2]
 
 ;define the L and Q parameters:
 Larr=2d9+(2.5d9-2d9)*dindgen(Nz)/(Nz-1)
 Qarr=0.01d0+(0.02d0-0.01d0)*dindgen(Nz)/(Nz-1)
            
 ;array of dimensions etc.:
 Lparms=[Nz,    $ ;number of voxels
         Nf,    $ ;number of frequencies
         NQ,    $ ;number of Q nodes
         NL,    $ ;number of L nodes
         NT,    $ ;number of temperatures
         0,     $ ;global DEM on/off key (on)
         0      ] ;global DDM on/off key (on)
 
 ;global floating-point parameters:
 Rparms=[1d18,  $ ;source area, cm^2
         3d9,   $ ;starting frequency, Hz 
         0.005  ] ;logarithmic step in frequency 

 ;single-voxel parameters:
 ParmLocal=dblarr(15)
 ParmLocal[0]=2d9    ;source depth, cm (total depth - the depths for individual voxels will be computed later)
 ParmLocal[1]=1d6    ;plasma temperature, K (not used in this example)
 ParmLocal[2]=1d9    ;electron/atomic concentration, cm^{-3} (not used in this example)
 ParmLocal[3]=100d0  ;magnetic field, G (will be changed later)
 ParmLocal[4]=120d0  ;viewing angle, degrees
 ParmLocal[5]=0d0    ;azimuthal angle, degrees
 ParmLocal[6]=0      ;emission mechanism flag (all on)
 ParmLocal[7]=30     ;maximum harmonic number
 ParmLocal[8]=0      ;proton concentration, cm^{-3} (not used in this example)
 ParmLocal[9]=0      ;neutral hydrogen concentration, cm^{-3}
 ParmLocal[10]=0     ;neutral helium concentration, cm^{-3}
 ParmLocal[11]=0     ;Q (will be changed later)
 ParmLocal[12]=0     ;L (will be changed later)
 ParmLocal[13]=0     ;element abundance code (coronal, following Feldman 1992)
 ParmLocal[14]=0     ;source area, cm^2 (not used in this example)
 
 ;multi-voxel parameters along line-of-sight:
 Parms=dblarr(15, Nz)
 for i=0, Nz-1 do begin
  Parms[*, i]=ParmLocal
  
  Parms[0, i]=ParmLocal[0]/Nz ;depths of individual voxels
  Parms[3, i]=1000d0-700d0*i/(Nz-1) ;magnetic field decreases linearly from 1000 to 300 G
  Parms[11, i]=Qarr[i]
  Parms[12, i]=Larr[i]
  
  ;all other parameters are the same in all voxels
 endfor
 
 RL=dblarr(7, Nf) ;output array
 
 ;main call:
 res=call_external(libname, 'GET_GX_MW', Lparms, Rparms, Parms, Qrun, Lrun, logtdem, DEM_cor_run, DDM_cor_run, RL)
 
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