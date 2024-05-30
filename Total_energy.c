#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Potential.h"

void Total_energy(double Total,double kinetic_energy,double R_ij,double alpha,double D,double Potential_energy,double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,double X[],double Y[],double Z[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],FILE *output1,int T,double d1,double d2,double d3,double boxs[],FILE *output8)
{
int i;
Potential(R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,output8,T);

Total=kinetic_energy+Potential_energy;

//printf("Total_energy %lf\n",Total);
fprintf(output1," %d   %lf\n",T,Total);

}
