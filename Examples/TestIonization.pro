pro TestIonization
 n0=1d13
 T0arr=exp(alog(1d3)+alog(1d5/1d3)*dindgen(100)/99)
 ne_arr=dblarr(100)
 nH_arr=dblarr(100)
 nHe_arr=dblarr(100)
 
 n_e=0d0
 n_H=0d0
 n_He=0d0
 
 libname='GRFF_DEM_Transfer_64.dll'
 
 for i=0, 99 do begin
  T0=T0arr[i]
  res=call_external(libname, 'FindIonizations', n0, T0, n_e, n_H, n_He, /d_value)
  ne_arr[i]=n_E
  nH_arr[i]=n_H
  nHe_arr[i]=n_He
 endfor
 res=call_external(libname, 'FindIonizations', n0, T0, n_e, n_H, n_He, /d_value, /unload)
 
 window, 4, title='Particle concentrations'
 wset, 4
 plot, T0arr, ne_arr, /xlog, /ylog, $
       xtitle='Temperature, K', ytitle='Concentration, cm!U-3!N'
 oplot, T0arr, nH_arr, linestyle=2
 oplot, T0arr, nHe_arr, linestyle=1
end