#include "MWtransfer.h"
#include "Coulomb.h"
#include "Zeta.h"
#include "Plasma.h"
#include "Messages.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv)
#else
extern "C" int GET_MW(int argc, void **argv)
#endif
{
 if (argc<6)
 {
  IDLmsg("GET_MW error: not enough parameters in the function call.");
  return -1;
 }

 int *Ndat=(int*)argv[0];
 double *Parms=(double*)argv[1];
 double *T_arr=(double*)argv[2];
 double *DEM_arr=(double*)argv[3];
 double *DDM_arr=(double*)argv[4];
 double *RL=(double*)argv[5];
           
 int res=MW_Transfer(Ndat[0], Ndat[1], Ndat[2], Ndat[3], Parms, T_arr, DEM_arr, DDM_arr, RL);

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) double CoulombLogarithm(int argc, void **argv)
#else
extern "C" double CoulombLogarithm(int argc, void **argv)
#endif
{
 double T=*((double*)argv[0]);
 double f=*((double*)argv[1]);

 return lnC1(T, f);
}

#ifndef LINUX
extern "C" __declspec(dllexport) double ZetaSolar(int argc, void **argv)
#else
extern "C" double ZetaSolar(int argc, void **argv)
#endif
{
 double T=*((double*)argv[0]);
 double f=*((double*)argv[1]);
 int ABcode=*((int*)argv[2]);

 return Zeta_Solar(T, f, ABcode);
}

#ifndef LINUX
extern "C" __declspec(dllexport) double SahaIonization(int argc, void **argv)
#else
extern "C" double SahaIonization(int argc, void **argv)
#endif
{
 double n0=*((double*)argv[0]);
 double T0=*((double*)argv[1]);
 return Saha(n0, T0);
}