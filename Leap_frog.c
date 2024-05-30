#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"force_calculation.h"

void Leap_frog(double alattice,int nx, int ny, int nz,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],int incremental_time,double X_Force[],double Y_Force[],double Z_Force[],int Totalatoms,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Dt,double mass,double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],int NSave,int T,double Vx_New[],double Vy_New[],double Vz_New[],double X_Old[],double Y_Old[],double Z_Old[],FILE *output4,double Old_Vel_X[],double Old_Vel_Y[],double Old_Vel_Z[],double New_Vel_X[],double New_Vel_Y[],double New_Vel_Z[],double d1,double d2,double d3,double boxs[],double Potential_energy)
{
int i;
double lx,ly,lz;
lx=alattice*nx;
ly=alattice*ny;
lz=alattice*nz;


fprintf(output4,"ITEM: TIMESTEP\n %d",T);
fprintf(output4,"\nITEM: NUMBER OF ATOMS\n%d",Totalatoms); 
fprintf(output4,"\nITEM: BOX BOUNDS");
fprintf(output4,"\n0.000 %lf xlo xhi\n",lx);
fprintf(output4,"0.000 %lf ylo yhi\n",ly);
fprintf(output4,"0.000 %lf zlo zhi\n",lz);
fprintf(output4,"\nITEM: ATOMS id x y z\n");

//previous position finding
for(i=0;i<Totalatoms;i++)
{
X_Old[i]=X[i]-(Vx_Scaled[i]*Dt)+(0.5*(Fx_ij[i]/mass)*Dt*Dt);

Y_Old[i]=Y[i]-(Vy_Scaled[i]*Dt)+(0.5*(Fy_ij[i]/mass)*Dt*Dt);

Z_Old[i]=Z[i]+(Vz_Scaled[i]*Dt)+(0.5*(Fz_ij[i]/mass)*Dt*Dt);

}

// updating forces
for(i=0;i<Totalatoms;i++)
{
X_Force[i]=Fx_ij[i];
Y_Force[i]=Fy_ij[i];
Z_Force[i]=Fz_ij[i];

}

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);

for(i=0;i<Totalatoms;i++)
{

Old_Vel_X[i]=(X[i]-X_Old[i])/Dt;
Old_Vel_Y[i]=(Y[i]-Y_Old[i])/Dt;
Old_Vel_Z[i]=(Z[i]-Z_Old[i])/Dt;


New_Vel_X[i]=Old_Vel_X[i]+(X_Force[i]*Dt)/mass;
New_Vel_Y[i]=Old_Vel_Y[i]+(Y_Force[i]*Dt)/mass;
New_Vel_Z[i]=Old_Vel_Z[i]+(Z_Force[i]*Dt)/mass;


X_New[i]=X[i]+(New_Vel_X[i]*Dt);
Y_New[i]=Y[i]+(New_Vel_Y[i]*Dt);
Z_New[i]=Z[i]+(New_Vel_Z[i]*Dt);

fprintf(output4,"%d %lf %lf %lf\n",i+1,X_New[i],Y_New[i],Z_New[i]);

Vx_Scaled[i]=(New_Vel_X[i]+Old_Vel_X[i])/2;
Vy_Scaled[i]=(New_Vel_Y[i]+Old_Vel_Y[i])/2;
Vz_Scaled[i]=(New_Vel_Z[i]+Old_Vel_Z[i])/2;

X_Old[i]=X[i];
Y_Old[i]=Y[i];
Z_Old[i]=Z[i];

X[i]=X_New[i];
Y[i]=Y_New[i];
Z[i]=Z_New[i];

}
}































