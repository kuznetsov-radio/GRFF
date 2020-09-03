#include <stdio.h>

const char* arr[]={
 "         dS;   1E+18   ;cm^2;                     Source/pixel area",
 "         dR;   1E+09   ;cm;                      Source/voxel depth",
 "        T_0;   1E+06   ;K;                       Plasma temperature",
 "        n_e;   1E+09   ;cm^{-3};                   Electron density",
 "          B;   200.0   ;G;                           Magnetic field",
 "      theta;    60.0   ;degrees;                      Viewing angle",
 "        psi;     0.0   ;degrees;                    Azimuthal angle",
 "      f_min;   1E+09   ;Hz;        Starting freq. to calc. spectrum",
 "         df;    0.02   ;Log(Hz);      Logarithmic step in frequency",
 "  mech_flag;       0   ;none;               Emission mechanism flag",
 "      s_max;      10   ;int;                Maximum harmonic number",
 "        n_H;       0   ;cm^{-3};     Neutral hydrogen concentration",
 "       n_He;       0   ;cm^{-3};       Neutral helium concentration",
 "    DEM_key;       0   ;0/1;                             DEM on/off",
 "  abund_key;       0   ;0/1/2;             Coronal / Caffau / Scott"
};

#define Nstrings 15

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS(int argc, void **argv)
#else
extern "C" float GET_PARMS(int argc, void **argv)
#endif
{
 FILE *F=fopen("Parms.txt", "w");
 if (F)
 {
  for (int i=0; i<Nstrings; i++) fprintf(F, "%s\n", arr[i]);
  fclose(F);
  return 0;
 }
 else return -1;
}