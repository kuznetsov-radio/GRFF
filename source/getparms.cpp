#include <stdio.h>

const char* arr1[]={
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "     N_temp;       0   ;int;     data;              Number of temperatures",
 "    DEM_key;       0   ;0/1;     user;         DEM on/off (for all voxels)",
 "    DDM_key;       0   ;0/1;     user;         DDM on/off (for all voxels)"
};

#define N1 5

const char* arr1_slice[]={
 "      N_pix;       1   ;int;     data;                    Number of pixels",
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "     N_temp;       0   ;int;     data;              Number of temperatures",
 "    DEM_key;       0   ;0/1;     user;         DEM on/off (for all voxels)",
 "    DDM_key;       0   ;0/1;     user;         DDM on/off (for all voxels)"
};

#define N1_slice 6

const char* arr1_extra_harmonics[]={
 "      s_min;       2   ;int;     user;             Minimum harmonic number",
 "      s_max;      10   ;int;     user;             Maximum harmonic number"
};

#define N1_extra_harmonics 2

const char* arr1_GX[]={
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "  N_Q_EBTEL;      60   ;int;     data;   No. of Q values in the EBTEL grid",
 "  N_L_EBTEL;      49   ;int;     data;   No. of L values in the EBTEL grid",
 "  N_T_EBTEL;     451   ;int;     data;       No. of temperatures (DEM/DDM)",
 "    DEM_key;       0   ;0/1;     user;         DEM on/off (for all voxels)",
 "    DDM_key;       0   ;0/1;     user;         DDM on/off (for all voxels)"
};

#define N1_GX 7

const char* arr1_GX_slice[]={
 "      N_pix;       1   ;int;     data;                    Number of pixels",
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "  N_Q_EBTEL;      60   ;int;     data;   No. of Q values in the EBTEL grid",
 "  N_L_EBTEL;      49   ;int;     data;   No. of L values in the EBTEL grid",
 "  N_T_EBTEL;     451   ;int;     data;       No. of temperatures (DEM/DDM)",
 "    DEM_key;       0   ;0/1;     user;         DEM on/off (for all voxels)",
 "    DDM_key;       0   ;0/1;     user;         DDM on/off (for all voxels)"
};

#define N1_GX_slice 8

const char* arr2[]={
 "         dS;   1E+18   ;cm^2;    data;                 Source/pixel area",
 "      f_min;   1E+09   ;Hz;      user;  Starting freq. to calc. spectrum",
 "         df;    0.02   ;Log(Hz); user;     Logarithmic step in frequency"
};

#define N2 3

const char* arr2_extra_harmonics[]={
 "    x_start;       0   ;cm;      data;                     LOS start (x)",
 "    y_start;       0   ;cm;      data;                     LOS start (y)",
 "    z_start;       0   ;cm;      data;                     LOS start (z)",
 "      x_end;       1   ;cm;      data;                       LOS end (x)",
 "      y_end;       1   ;cm;      data;                       LOS end (y)",
 "      z_end;       1   ;cm;      data;                       LOS end (z)"
};

#define N2_extra_harmonics 6

const char* arr3[]={
 "         dR;   1E+09   ;cm;      data;                Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;     Electron or total gas density",
 "          B;   200.0   ;G;       data;                    Magnetic field",
 "      theta;    60.0   ;degrees; data;                     Viewing angle",
 "        phi;     0.0   ;degrees; data;                   Azimuthal angle",
 "  mech_flag;       0   ;none;    data;           Emission mechanism flag",
 "      s_max;      10   ;int;     data;           Maximum harmonic number",
 "        n_p;       0   ;cm^{-3}; data;              Proton concentration",
 "       n_HI;       0   ;cm^{-3}; data;    Neutral hydrogen concentration",
 "      n_HeI;       0   ;cm^{-3}; data;      Neutral helium concentration",
 "DEM_key_loc;       0   ;0/1;     data;                        DEM on/off",
 "DDM_key_loc;       0   ;0/1;     data;                        DDM on/off",
 "  abund_key;       0   ;0/1/2;   data;          Coronal / Caffau / Scott",
 "      S_loc;       0   ;cm^2;    data;                 Local source area"
};

#define N3 15

const char* arr3_harmonics[]={
 "         dR;   1E+09   ;cm;      data;                Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;     Electron or total gas density",
 "          B;   200.0   ;G;       data;                    Magnetic field",
 "      theta;    60.0   ;degrees; data;                     Viewing angle",
 "        phi;     0.0   ;degrees; data;                   Azimuthal angle",
 "  mech_flag;       0   ;none;    data;           Emission mechanism flag",
 "     unused;      10   ;int;     data;                            unused",
 "        n_p;       0   ;cm^{-3}; data;              Proton concentration",
 "       n_HI;       0   ;cm^{-3}; data;    Neutral hydrogen concentration",
 "      n_HeI;       0   ;cm^{-3}; data;      Neutral helium concentration",
 "DEM_key_loc;       0   ;0/1;     data;                        DEM on/off",
 "DDM_key_loc;       0   ;0/1;     data;                        DDM on/off",
 "  abund_key;       0   ;0/1/2;   data;          Coronal / Caffau / Scott",
 "      S_loc;       0   ;cm^2;    data;                 Local source area"
};

#define N3_harmonics 15

const char* arr3_GX[]={
 "         dR;   1E+09   ;cm;      data;                Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;     Electron or total gas density",
 "          B;   200.0   ;G;       data;                    Magnetic field",
 "      theta;    60.0   ;degrees; data;                     Viewing angle",
 "        phi;     0.0   ;degrees; data;                   Azimuthal angle",
 "  mech_flag;       0   ;none;    data;           Emission mechanism flag",
 "      s_max;      10   ;int;     data;           Maximum harmonic number",
 "        n_p;       0   ;cm^{-3}; data;              Proton concentration",
 "       n_HI;       0   ;cm^{-3}; data;    Neutral hydrogen concentration",
 "      n_HeI;       0   ;cm^{-3}; data;      Neutral helium concentration",
 "          Q;       0   ;relative;data;                      Heating rate",
 "          L;       0   ;cm;      data;                       Loop length",
 "  abund_key;       0   ;0/1/2;   data;          Coronal / Caffau / Scott",
 "      S_loc;       0   ;cm^2;    data;                 Local source area"
};

#define N3_GX 15

void WriteParms(const char **arr, const char *fname, int N, int add)
{
 FILE *F=fopen(fname, add ? "a" : "w");
 if (F)
 {
  for (int i=0; i<N; i++) fprintf(F, "%s\n", arr[i]);
  fclose(F);
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS(int argc, void **argv) //standard single-thread
#else
extern "C" float GET_PARMS(int argc, void **argv)
#endif
{
 WriteParms(arr1, "Long_input.txt", N1, 0);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr3, "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS_SLICE(int argc, void **argv) //standard multi-thread
#else
extern "C" float GET_PARMS_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1_slice, "Long_input.txt", N1_slice, 0);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr3, "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS1(int argc, void **argv) //gyrolayers single-thread
#else
extern "C" float GET_PARMS1(int argc, void **argv)
#endif
{
 WriteParms(arr1, "Long_input.txt", N1, 0);
 WriteParms(arr1_extra_harmonics, "Long_input.txt", N1_extra_harmonics, 1);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr2_extra_harmonics, "Real_input.txt", N2_extra_harmonics, 1);
 WriteParms(arr3_harmonics, "Parms_input.txt", N3_harmonics, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS1_SLICE(int argc, void **argv) //gyrolayers multi-thread
#else
extern "C" float GET_PARMS1_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1_slice, "Long_input.txt", N1_slice, 0);
 WriteParms(arr1_extra_harmonics, "Long_input.txt", N1_extra_harmonics, 1);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr2_extra_harmonics, "Real_input.txt", N2_extra_harmonics, 1);
 WriteParms(arr3_harmonics, "Parms_input.txt", N3_harmonics, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_GX_PARMS(int argc, void **argv) //EBTEL-based single-thread
#else
extern "C" float GET_GX_PARMS(int argc, void **argv)
#endif
{
 WriteParms(arr1_GX, "Long_input.txt", N1_GX, 0);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr3_GX, "Parms_input.txt", N3_GX, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_GX_PARMS_SLICE(int argc, void **argv) //EBTEL-based multi-thread
#else
extern "C" float GET_GX_PARMS_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1_GX_slice, "Long_input.txt", N1_GX_slice, 0);
 WriteParms(arr2, "Real_input.txt", N2, 0);
 WriteParms(arr3_GX, "Parms_input.txt", N3_GX, 0);
 return 0;
}