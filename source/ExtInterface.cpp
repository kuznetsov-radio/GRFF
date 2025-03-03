#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "MWtransfer.h"
#include "Coulomb.h"
#include "Zeta.h"
#include "Plasma.h"
#include "ExtInterface.h"
#include "GXcode.h"
#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_main(int argc, void **argv)
#else
extern "C" int GET_MW_main(int argc, void **argv)
#endif
{
 int res=0;

 if (argc==7 || argc==10)
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms=(double*)argv[2];
  double *T_arr=(double*)argv[3];
  double *DEM_arr=(double*)argv[4];
  double *DDM_arr=(double*)argv[5];

  if (argc==7)
  {
   double *RL=(double*)argv[6];
           
   res=MW_Transfer(Lparms, Rparms, Parms, T_arr, DEM_arr, DDM_arr, RL, 0, 0, 0, 0, 0, 0);
  }
  else
  {
   double *fZ_arr=(double*)argv[6];
   double *TZ_arr=(double*)argv[7];
   double *Z_arr=(double*)argv[8];
   double *RL=(double*)argv[9];

   res=MW_Transfer(Lparms, Rparms, Parms, T_arr, DEM_arr, DDM_arr, RL, 1, fZ_arr, TZ_arr, Z_arr, 0, 0);
  }
 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv)
#else
extern "C" int GET_MW(int argc, void **argv)
#endif
{
 int res=0;

 if (argc==7 || argc==10)
 {
  int *Lparms1=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms1=(double*)argv[2];
  double *T_arr=(double*)argv[3];
  double *DEM_arr=(double*)argv[4];
  double *DDM_arr=(double*)argv[5];
  
  #define InSize1 15

  int Nz=Lparms1[0];
  double *Parms=(double*)malloc(sizeof(double)*InSize*Nz);

  int NT=Lparms1[2];

  int DEM_on_global=(Lparms1[3]==0);
  int DDM_on_global=(Lparms1[4]==0);

  int Lparms[6];
  for (int i=0; i<3; i++) Lparms[i]=Lparms1[i];
  if (argc==10) for (int i=3; i<6; i++) Lparms[i]=Lparms1[i+2];

  for (int j=0; j<Nz; j++)
  {
   double *p=Parms+j*InSize;
   double *p1=Parms1+j*InSize1;

   for (int i=0; i<=7; i++) p[i]=p1[i]; //parameters 0-7 are the same in both functions
   p[8]=p1[9]; //n_H
   p[9]=p1[10]; //n_He
   p[10]=DEM_on_global ? p1[11] : 1; //DEM key
   p[11]=DDM_on_global ? p1[12] : 1; //DDM key
   p[12]=p1[13]; //abundance key
   p[13]=0; //Maxwellian distribution only
   p[14]=p1[14]; //source area (if nonzero)

   int DEM_on=(p[10]==0 && NT>1);
   int DDM_on=(p[11]==0 && NT>1);

   if (!DEM_on && !DDM_on)
   {
    double T0=p1[1];

    if (T0<1e5)
    {
     double n_p=p1[8];
     double n_H=p1[9];

     if (n_p==0 && n_H==0)
     {
      double n0=p1[2];

      double n_e, n_He;

      FindIonizationsSolar(n0, T0, &n_e, &n_H, &n_He);

      p[2]=n_e;
      p[8]=n_H;
      p[9]=n_He; 
     }
    } //otherwise, if T0=>1e5, n_e (p[2]) has already been assigned equal to n_0(p1[2])
   }
  }

  if (argc==7)
  {
   double *RL=(double*)argv[6];

   res=MW_Transfer(Lparms, Rparms, Parms, T_arr, DEM_arr, DDM_arr, RL, 0, 0, 0, 0, 0, 0);
  }
  else
  {
   double *fZ_arr=(double*)argv[6];
   double *TZ_arr=(double*)argv[7];
   double *Z_arr=(double*)argv[8];
   double *RL=(double*)argv[9];

   res=MW_Transfer(Lparms, Rparms, Parms, T_arr, DEM_arr, DDM_arr, RL, 1, fZ_arr, TZ_arr, Z_arr, 0, 0);
  }

  free(Parms);

 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_SLICE(int argc, void** argv)
#else
extern "C" int GET_MW_SLICE(int argc, void** argv)
#endif
{
 if (argc<7) return -1;

 int *Lparms_M=(int*)argv[0];
 double *Rparms_M=(double*)argv[1];
 double *Parms_M=(double*)argv[2];
 double *T_arr=(double*)argv[3];
 double *DEM_arr_M=(double*)argv[4];
 double *DDM_arr_M=(double*)argv[5];
 double *RL_M=(double*)argv[6];

 int Npix=Lparms_M[0];
 int Nz=Lparms_M[1];
 int Nf=Lparms_M[2];
 int NT=Lparms_M[3];

 #ifndef LINUX

 concurrency::parallel_for(0, Npix, [&](int pix)
 {
  void *ARGV[7];
  ARGV[0]=(void*)(Lparms_M+1);
  ARGV[1]=(void*)(Rparms_M+pix*RpSize);
  ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
  ARGV[3]=(void*)T_arr;
  ARGV[4]=(void*)(DEM_arr_M+pix*Nz*NT);
  ARGV[5]=(void*)(DDM_arr_M+pix*Nz*NT);
  ARGV[6]=(void*)(RL_M+pix*Nf*OutSize);

  GET_MW(7, ARGV);
 });

 #else

 #pragma omp parallel for
 for(int pix=0; pix<Npix; pix++)
 {
  void *ARGV[7];
  ARGV[0]=(void*)(Lparms_M+1);
  ARGV[1]=(void*)(Rparms_M+pix*RpSize);
  ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
  ARGV[3]=(void*)T_arr;
  ARGV[4]=(void*)(DEM_arr_M+pix*Nz*NT);
  ARGV[5]=(void*)(DDM_arr_M+pix*Nz*NT);
  ARGV[6]=(void*)(RL_M+pix*Nf*OutSize);

  GET_MW(7, ARGV);
 }

 #endif

 return 0;
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
extern "C" __declspec(dllexport) double SahaIonizationH(int argc, void **argv)
#else
extern "C" double SahaIonizationH(int argc, void **argv)
#endif
{
 double n0=*((double*)argv[0]);
 double T0=*((double*)argv[1]);
 return SahaH(n0, T0);
}

#ifndef LINUX
extern "C" __declspec(dllexport) double FindIonizations(int argc, void **argv)
#else
extern "C" double FindIonizations(int argc, void **argv)
#endif
{
 double n0=*((double*)argv[0]);
 double T0=*((double*)argv[1]);
 double *n_e=(double*)argv[2];
 double *n_H=(double*)argv[3];
 double *n_He=(double*)argv[4];

 FindIonizationsSolar(n0, T0, n_e, n_H, n_He);

 return *n_e/n0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW1(int argc, void **argv)
#else
extern "C" int GET_MW1(int argc, void **argv)
#endif
{
 if (argc!=8) return -1;
 else
 { 
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms1=(double*)argv[2];
  double *T_arr=(double*)argv[3];
  double *DEM_arr=(double*)argv[4];
  double *DDM_arr=(double*)argv[5];
  double *RL=(double*)argv[6];
  double *GRparms=(double*)argv[7];
  
  int *srange=Lparms+5;
    
  #define InSize1 15

  int Nz=Lparms[0];
  double *Parms=(double*)malloc(sizeof(double)*InSize*Nz);

  int NT=Lparms[2];

  int DEM_on_global=(Lparms[3]==0);
  int DDM_on_global=(Lparms[4]==0);

  for (int j=0; j<Nz; j++)
  {
   double *p=Parms+j*InSize;
   double *p1=Parms1+j*InSize1;

   for (int i=0; i<=7; i++) p[i]=p1[i]; //parameters 0-7 are the same in both functions
   p[8]=p1[9]; //n_H
   p[9]=p1[10]; //n_He
   p[10]=DEM_on_global ? p1[11] : 1; //DEM key
   p[11]=DDM_on_global ? p1[12] : 1; //DDM key
   p[12]=p1[13]; //abundance key
   p[13]=0; //Maxwellian distribution only
   p[14]=p1[14]; //source area (if nonzero)

   int DEM_on=(p[10]==0 && NT>1);
   int DDM_on=(p[11]==0 && NT>1);

   if (!DEM_on && !DDM_on)
   {
    double T0=p1[1];

    if (T0<1e5)
    {
     double n_p=p1[8];
     double n_H=p1[9];

     if (n_p==0 && n_H==0)
     {
      double n0=p1[2];

      double n_e, n_He;

      FindIonizationsSolar(n0, T0, &n_e, &n_H, &n_He);

      p[2]=n_e;
      p[8]=n_H;
      p[9]=n_He; 
     }
    } //otherwise, if T0=>1e5, n_e (p[2]) has already been assigned equal to n_0(p1[2])
   }
  }

  int res=MW_Transfer(Lparms, Rparms, Parms, T_arr, DEM_arr, DDM_arr, RL, 0, 0, 0, 0, srange, GRparms);

  free(Parms);

  return res;
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW1_SLICE(int argc, void** argv)
#else
extern "C" int GET_MW1_SLICE(int argc, void** argv)
#endif
{
 if (argc!=8) return -1;
 else
 {
  int *Lparms_M=(int*)argv[0];
  double *Rparms_M=(double*)argv[1];
  double *Parms_M=(double*)argv[2];
  double *T_arr=(double*)argv[3];
  double *DEM_arr_M=(double*)argv[4];
  double *DDM_arr_M=(double*)argv[5];
  double *RL_M=(double*)argv[6];
  double *GRparms_M=(double*)argv[7];

  int *srange=Lparms_M+6;
  
  int Npix=Lparms_M[0];
  int Nz=Lparms_M[1];
  int Nf=Lparms_M[2];
  int NT=Lparms_M[3];
  int Ns=srange ? srange[1]-srange[0]+1 : 0;
  if (Ns<0) Ns=0;

  #ifndef LINUX

  concurrency::parallel_for(0, Npix, [&](int pix)
  {
   void *ARGV[8];
   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSize1);
   ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
   ARGV[3]=(void*)T_arr;
   ARGV[4]=(void*)(DEM_arr_M+pix*Nz*NT);
   ARGV[5]=(void*)(DDM_arr_M+pix*Nz*NT);
   ARGV[6]=(void*)(RL_M+pix*Nf*OutSize);
   ARGV[7]=Ns ? (void*)(GRparms_M+pix*GpSize*Nf*Ns) : 0;

   GET_MW1(8, ARGV);
  });

  #else

  #pragma omp parallel for
  for(int pix=0; pix<Npix; pix++)
  {
   void *ARGV[8];
   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSize1);
   ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
   ARGV[3]=(void*)T_arr;
   ARGV[4]=(void*)(DEM_arr_M+pix*Nz*NT);
   ARGV[5]=(void*)(DDM_arr_M+pix*Nz*NT);
   ARGV[6]=(void*)(RL_M+pix*Nf*OutSize);
   ARGV[7]=Ns ? (void*)(GRparms_M+pix*GpSize*Nf*Ns) : 0;

   GET_MW1(8, ARGV);
  }

  #endif

  return 0;
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_GX_MW(int argc, void **argv)
#else
extern "C" int GET_GX_MW(int argc, void **argv)
#endif
{
 int res=0;

 if (argc>=9)
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms=(double*)argv[2];
  float *Qrun=(float*)argv[3];
  float *Lrun=(float*)argv[4];
  float *logT=(float*)argv[5];
  float *DEMrun=(float*)argv[6];
  float *DDMrun=(float*)argv[7];
  double *RL=(double*)argv[8];

  int Nz=Lparms[0];
  int NQ=Lparms[2];
  int NL=Lparms[3];
  int NT=Lparms[4];

  int Lparms_o[LpSize]={Lparms[0], Lparms[1], Lparms[4], Lparms[5], Lparms[6]};

  double *T_arr=(logT!=0 && NT>0) ? (double*)malloc(sizeof(double)*NT) : 0;
  if (T_arr) for (int i=0; i<NT; i++) T_arr[i]=pow(10.0, double(logT[i]));

  double *Parms_o=(double*)malloc(sizeof(double)*Nz*InSize);
  memcpy(Parms_o, Parms, sizeof(double)*Nz*InSize);

  double *DEM_arr=(DEMrun!=0 && Qrun!=0 && Lrun!=0 && NT>0 && NQ>0 && NL>0) ? (double*)malloc(sizeof(double)*NT*Nz) : 0;
  double *DDM_arr=(DDMrun!=0 && Qrun!=0 && Lrun!=0 && NT>0 && NQ>0 && NL>0) ? (double*)malloc(sizeof(double)*NT*Nz) : 0;

  for (int i=0; i<Nz; i++)
  {
   double Q=Parms[D2(InSize, 11, i)];
   double L=Parms[D2(InSize, 12, i)];

   if (Q>0 && L>0 && (DEM_arr!=0 || DDM_arr!=0))
   {
    int ok=InterpolateEBTEL(NQ, NL, NT, Q, L, Qrun, Lrun, DEMrun, DDMrun, 
                            DEM_arr ? DEM_arr+i*NT : 0, DDM_arr ? DDM_arr+i*NT : 0);

    Parms_o[D2(InSize, 11, i)]=!(DEM_arr!=0 && ok!=0);
    Parms_o[D2(InSize, 12, i)]=!(DDM_arr!=0 && ok!=0);
   }
   else Parms_o[D2(InSize, 11, i)]=Parms_o[D2(InSize, 12, i)]=1;
  }

  void *ARGV[7];
  ARGV[0]=(void*)Lparms_o;
  ARGV[1]=(void*)Rparms;
  ARGV[2]=(void*)Parms_o;
  ARGV[3]=(void*)T_arr;
  ARGV[4]=(void*)DEM_arr;
  ARGV[5]=(void*)DDM_arr;
  ARGV[6]=(void*)RL;

  res=GET_MW(7, ARGV);

  if (DDM_arr) free(DDM_arr);
  if (DEM_arr) free(DEM_arr);
  free(Parms_o);
  if (T_arr) free(T_arr);
 }
 else res=-1;

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_GX_MW_SLICE(int argc, void **argv)
#else
extern "C" int GET_GX_MW_SLICE(int argc, void **argv)
#endif
{
 if (argc>=9)
 {
  int *Lparms_M=(int*)argv[0];
  double *Rparms_M=(double*)argv[1];
  double *Parms_M=(double*)argv[2];
  float *Qrun=(float*)argv[3];
  float *Lrun=(float*)argv[4];
  float *logT=(float*)argv[5];
  float *DEMrun=(float*)argv[6];
  float *DDMrun=(float*)argv[7];
  double *RL_M=(double*)argv[8];

  int Npix=Lparms_M[0];
  int Nz=Lparms_M[1];
  int Nf=Lparms_M[2];

  #ifndef LINUX
  concurrency::parallel_for(0, Npix, [&](int pix)
  #else
  #pragma omp parallel for
  for(int pix=0; pix<Npix; pix++)
  #endif
  {
  void *ARGV[9];
  ARGV[0]=(void*)(Lparms_M+1);
  ARGV[1]=(void*)(Rparms_M+pix*RpSize);
  ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
  ARGV[3]=(void*)Qrun;
  ARGV[4]=(void*)Lrun;
  ARGV[5]=(void*)logT;
  ARGV[6]=(void*)DEMrun;
  ARGV[7]=(void*)DDMrun;
  ARGV[8]=(void*)(RL_M+pix*Nf*OutSize);

  GET_GX_MW(9, ARGV);
  #ifndef LINUX
  });
  #else
  }
  #endif

  return 0;
 }
 else return -1;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int PyGET_MW_USER(int *Lparms, double *Rparms, double *Parms,
                                                   double *T_arr,double *DEM_arr, double *DDM_arr, 
                                                   double *fzeta_arr, double *Tzeta_arr, double *zeta_arr,
                                                   double *RL)
#else
extern "C" int PyGET_MW_USER(int *Lparms, double *Rparms, double *Parms,
                           double *T_arr,double *DEM_arr, double *DDM_arr, 
                           double *fzeta_arr, double *Tzeta_arr, double *zeta_arr,
                           double *RL)
#endif
{
 void *ARGV[10];
 ARGV[0]=(void*)Lparms;
 ARGV[1]=(void*)Rparms;
 ARGV[2]=(void*)Parms;
 ARGV[3]=(void*)T_arr;
 ARGV[4]=(void*)DEM_arr;
 ARGV[5]=(void*)DDM_arr;
 ARGV[6]=(void*)fzeta_arr;
 ARGV[7]=(void*)Tzeta_arr;
 ARGV[8]=(void*)zeta_arr;
 ARGV[9]=(void*)RL; 

 return GET_MW(10, ARGV);
}


#ifndef LINUX
extern "C" __declspec(dllexport) int PyGET_MW(int *Lparms, double *Rparms, double *Parms,
                                              double *T_arr,double *DEM_arr, double *DDM_arr,
                                              double *RL)
#else
extern "C" int PyGET_MW(int *Lparms, double *Rparms, double *Parms,
                        double *T_arr,double *DEM_arr, double *DDM_arr,
                        double *RL)
#endif
{
 void *ARGV[7];
 ARGV[0]=(void*)Lparms;
 ARGV[1]=(void*)Rparms;
 ARGV[2]=(void*)Parms;
 ARGV[3]=(void*)T_arr;
 ARGV[4]=(void*)DEM_arr;
 ARGV[5]=(void*)DDM_arr;
 ARGV[6]=(void*)RL; 

 return GET_MW(7, ARGV);
}
