#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"force_calculation.h"

void Verlet(double alattice,int nx, int ny, int nz,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],int incremental_time,double X_Force[],double Y_Force[],double Z_Force[],int Totalatoms,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Dt,double mass,double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],FILE *output,int NSave,int T,double Vx_New[],double Vy_New[],double Vz_New[],double X_Old[],double Y_Old[],double Z_Old[],FILE *output3,double d1,double d2,double d3,double boxs[],double Potential_energy)
{
int i;
double lx,ly,lz;
lx=alattice*nx;
ly=alattice*ny;
lz=alattice*nz;

printf("mass is %lf \n",mass);
printf("Dt is %lf \n",Dt);
//if(T%NSave==0)
//{

fprintf(output3,"ITEM: TIMESTEP\n %d",T);
//fprintf(output," %d\n",incremental_time);
fprintf(output3,"\nITEM: NUMBER OF ATOMS\n%d",Totalatoms); 
fprintf(output3,"\nITEM: BOX BOUNDS");
fprintf(output3,"\n0.000 %lf xlo xhi\n",lx);
fprintf(output3,"0.000 %lf ylo yhi\n",ly);
fprintf(output3,"0.000 %lf zlo zhi\n",lz);
fprintf(output3,"\nITEM: ATOMS id x y z\n");



//previous position finding
for(i=0;i<Totalatoms;i++)
{
X_Old[i]=X[i]-(Vx_Scaled[i]*Dt)+(0.5*(Fx_ij[i]/mass)*Dt*Dt);

Y_Old[i]=Y[i]-(Vy_Scaled[i]*Dt)+(0.5*(Fy_ij[i]/mass)*Dt*Dt);

Z_Old[i]=Z[i]+(Vz_Scaled[i]*Dt)+(0.5*(Fz_ij[i]/mass)*Dt*Dt);



//fprintf(output,"%d %lf %lf %lf\n",i+1,X_New[i],Y_New[i],Z_New[i]);
}

//updating position
for(i=0;i<Totalatoms;i++)
{
X_New[i]= 2*X[i]-X_Old[i];
X[i]=X_New[i];
Y_New[i]= 2*Y[i]-Y_Old[i];
Y[i]=Y_New[i];
Z_New[i]= 2*Z[i]-Z_Old[i];
Z[i]=Z_New[i];

fprintf(output3,"%d %lf %lf %lf\n",i+1,X_New[i],Y_New[i],Z_New[i]);
}

//Velocity update
for(i=0;i<Totalatoms;i++)
{
Vx_Scaled[i]=(X_New[i]-X_Old[i])/(2*Dt);
Vy_Scaled[i]=(Y_New[i]-Y_Old[i])/(2*Dt);
Vz_Scaled[i]=(Z_New[i]-Z_Old[i])/(2*Dt);
}


// updating forces
for(i=0;i<Totalatoms;i++)
{
X_Force[i]=Fx_ij[i];
Y_Force[i]=Fy_ij[i];
Z_Force[i]=Fz_ij[i];

}

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);




//}

}
