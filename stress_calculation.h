#ifndef STRESS_CALCULATION_H_INCLUDED
#define STRESS_CALCULATION_H_INCLUDED

void stress_calculation(double d1,double d2,double d3,double Fx_ij[],double Fy_ij[],double Fz_ij[],double Force_y,double stress,FILE *output5,int T,int Totalatoms,double boxs[],int s);
#endif
