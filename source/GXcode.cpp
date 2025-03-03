#include <math.h>
#include <memory.h>
#include <float.h>
#include "ExtMath.h"
#include "ExtInterface.h"

int InterpolateEBTEL(int NQ, int NL, int NT, double Q, double L, 
	                 float *Qrun, float *Lrun, float *DEM_run, float *DDM_run, 
	                 double *DEM, double *DDM)
{
 int res=0;

 if (DEM!=0) memset(DEM, 0, sizeof(double)*NT);
 if (DDM!=0) memset(DDM, 0, sizeof(double)*NT);

 int Lind=value_locate_M(Lrun, NL, NQ, L);

 if (Lind>=0 && Lind<(NL-1))
 {
  int Qind1=value_locate_M(Qrun+Lind*NQ, NQ, 1, Q);
  int Qind2=value_locate_M(Qrun+(Lind+1)*NQ, NQ, 1, Q);

  if (Qind1>=0 && Qind1<(NQ-1) && Qind2>=0 && Qind2<(NQ-1))
  {
   res=1;

   double Llog=log(L);
   double Lgrid1=log(double(Lrun[D2(NQ, 0, Lind)]));
   double Lgrid2=log(double(Lrun[D2(NQ, 0, Lind+1)]));
   double dL=(Llog-Lgrid1)/(Lgrid2-Lgrid1);

   double Qlog=log(Q);
   double Qgrid11=log(double(Qrun[D2(NQ, Qind1, Lind)]));
   double Qgrid12=log(double(Qrun[D2(NQ, Qind1+1, Lind)]));
   double Qgrid21=log(double(Qrun[D2(NQ, Qind2, Lind+1)]));
   double Qgrid22=log(double(Qrun[D2(NQ, Qind2+1, Lind+1)]));
   double dQ1=(Qlog-Qgrid11)/(Qgrid12-Qgrid11);
   double dQ2=(Qlog-Qgrid21)/(Qgrid22-Qgrid21);

   if (DEM!=0) for (int l=0; l<NT; l++)
   {
    DEM[l]=DEM_run[D3(NT, NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
           DEM_run[D3(NT, NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
           DEM_run[D3(NT, NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
           DEM_run[D3(NT, NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
	if (!finite(DEM[l])) DEM[l]=0;
   }

   if (DDM!=0) for (int l=0; l<NT; l++)
   {
    DDM[l]=DDM_run[D3(NT, NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
           DDM_run[D3(NT, NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
           DDM_run[D3(NT, NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
           DDM_run[D3(NT, NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
	if (!finite(DDM[l])) DDM[l]=0;
   }
  }
 }

 return res;
}