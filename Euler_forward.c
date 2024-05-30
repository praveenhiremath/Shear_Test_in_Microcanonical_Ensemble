#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"force_calculation.h"

void Euler_forward(FILE *output2,double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,int incremental_time,double mass,double Dt,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],double Vx_New[],double Vy_New[],double Vz_New[],int T,double d1,double d2,double d3,double boxs[],double Potential_energy)
{
int i;
double lx,ly,lz;
lx=alattice*nx;
ly=alattice*ny;
lz=alattice*nz;


fprintf(output2,"ITEM: TIMESTEP\n %d",T);
fprintf(output2,"\nITEM: NUMBER OF ATOMS\n%d",Totalatoms); 
fprintf(output2,"\nITEM: BOX BOUNDS");
fprintf(output2,"\n0.000 %lf xlo xhi\n",lx);
fprintf(output2,"0.000 %lf ylo yhi\n",ly);
fprintf(output2,"0.000 %lf zlo zhi\n",lz);
fprintf(output2,"\nITEM: ATOMS id x y z\n");



//updating position
for(i=0;i<Totalatoms;i++)
{
X_New[i]=X[i]+(Dt*Vx_New[i])+(0.5*(Fx_ij[i]/mass)*Dt*Dt);
X[i]=X_New[i];
Y_New[i]=Y[i]+(Dt*Vy_New[i])+(0.5*(Fy_ij[i]/mass)*Dt*Dt);
Y[i]=Y_New[i];
Z_New[i]=Z[i]+(Dt*Vz_New[i])+(0.5*(Fz_ij[i]/mass)*Dt*Dt);
Z[i]=Z_New[i];
fprintf(output2,"%d %lf %lf %lf\n",i+1,X_New[i],Y_New[i],Z_New[i]);
}

//Velocity updating

for(i=0;i<Totalatoms;i++)
{

Vx_Scaled[i]=Vx_Scaled[i]+((Dt)*(Fx_ij[i]/mass));
Vx_New[i]=Vx_Scaled[i];
Vy_Scaled[i]=Vy_Scaled[i]+((Dt)*(Fy_ij[i]/mass));
Vy_New[i]=Vy_Scaled[i];
Vz_Scaled[i]=Vz_Scaled[i]+((Dt)*(Fz_ij[i]/mass));
Vz_New[i]=Vz_Scaled[i];
//printf("\n Velocities updated %lf %lf %lf\n",Vx_Scaled[i],Vy_Scaled[i],Vz_Scaled[i]);
//printf("Vy_Scaled %lf\n",Vy_Scaled[i]);
}

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);



}


