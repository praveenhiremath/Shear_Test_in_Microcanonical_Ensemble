#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"force_calculation.h"

void Potential(double R_ij,double alpha,double D,double Potential_energy,double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,double X[],double Y[],double Z[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],double d1,double d2,double d3,double boxs[],FILE *output8,int T)
{
double A1,A2,R_0;
R_0=3.4068;
force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);
A1=exp((-2*alpha)*(R_ij-R_0));
A2=exp(((-alpha)*(R_ij-R_0)));
Potential_energy= (D)*(A1-(2*A2));
//printf("Potential is %lf\n",Potential_energy);
fprintf(output8,"\n %d\t %lf\n",T,Potential_energy);

}
