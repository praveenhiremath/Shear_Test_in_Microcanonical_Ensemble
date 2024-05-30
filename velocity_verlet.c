#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"force_calculation.h"

void velocity_verlet(double alattice,int nx, int ny, int nz,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],int incremental_time,double X_Force[],double Y_Force[],double Z_Force[],int Totalatoms,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Dt,double mass,double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],FILE *output,int NSave,int T,double Vx_New[],double Vy_New[],double Vz_New[],double V_scaled,double d1,double d2,double d3,double boxs[],double Potential_energy)
{
int i;
double lx,ly,lz;
lx=alattice*nx;
ly=alattice*ny;
lz=alattice*nz;


//if(T%NSave==0)
//{

fprintf(output,"ITEM: TIMESTEP\n %d",T);
//fprintf(output," %d\n",incremental_time);
fprintf(output,"\nITEM: NUMBER OF ATOMS\n%d",Totalatoms); 
fprintf(output,"\nITEM: BOX BOUNDS");
fprintf(output,"\n0.000 %lf xlo xhi\n",lx);
fprintf(output,"0.000 %lf ylo yhi\n",ly);
fprintf(output,"0.000 %lf zlo zhi\n",lz);
fprintf(output,"\nITEM: ATOMS id x y z\n");



//updating position
for(i=0;i<Totalatoms;i++)
{
X_New[i]=X[i]+(Vx_New[i]*Dt)+(0.5*(Fx_ij[i]/mass)*Dt*Dt);
X[i]=X_New[i];
Y_New[i]=Y[i]+(Vy_New[i]*Dt)+(0.5*(Fy_ij[i]/mass)*Dt*Dt);
Y[i]=Y_New[i];
Z_New[i]=Z[i]+(Vz_New[i]*Dt)+(0.5*(Fz_ij[i]/mass)*Dt*Dt);
Z[i]=Z_New[i];

//X[i]=X[i]+(Vx_Scaled[i]*Dt)+(0.5*(Fx_ij[i]/mass)*Dt*Dt);

//Y[i]=Y[i]+(Vx_Scaled[i]*Dt)+(0.5*(Fy_ij[i]/mass)*Dt*Dt);

//Z[i]=Z[i]+(Vx_Scaled[i]*Dt)+(0.5*(Fz_ij[i]/mass)*Dt*Dt);

fprintf(output,"%d %lf %lf %lf\n",i+1,X_New[i],Y_New[i],Z_New[i]);
}
// updating forces
for(i=0;i<Totalatoms;i++)
{
X_Force[i]=Fx_ij[i];
Y_Force[i]=Fy_ij[i];
Z_Force[i]=Fz_ij[i];

}


force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);
//Updating velocities
for(i=0;i<Totalatoms;i++)
{
Vx_Scaled[i]=Vx_Scaled[i]+((0.5*Dt)*((Fx_ij[i]+X_Force[i])/mass));
Vx_New[i]=Vx_Scaled[i];
Vy_Scaled[i]=Vy_Scaled[i]+((0.5*Dt)*((Fy_ij[i]+Y_Force[i])/mass));
Vy_New[i]=Vy_Scaled[i];
Vz_Scaled[i]=Vz_Scaled[i]+((0.5*Dt)*((Fz_ij[i]+Z_Force[i])/mass));
Vz_New[i]=Vz_Scaled[i];
//printf("\n Velocities updated %lf %lf %lf\n",Vx_Scaled[i],Vy_Scaled[i],Vz_Scaled[i]);
}

V_scaled=0.0;
for(i = 0 ; i <Totalatoms ; i++)
		    {
		        //printf("%d %lf %lf %lf\n",i+1, Vx_Scaled[i], Vy_Scaled[i], Vz_Scaled[i]);
//Sum of sqaured velocities of atoms                         
                        V_scaled=V_scaled+(Vx_Scaled[i]*Vx_Scaled[i])+(Vy_Scaled[i]*Vy_Scaled[i])+(Vz_Scaled[i]*Vz_Scaled[i]);
                        
//                      printf("Scaled Velocity is %lf\n",V_scaled);
		    }

//printf("Scaled Velocity is %lf\n",V_scaled);

//}

}
