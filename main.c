//Praveenkumar Hiremath 
//Matriculation Number: 57955
//Assignment 3: Modelling a shear test in the micro canonical ensemble 


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"bcc_structure.h"
#include"fcc_structure.h"
#include"hcp_structure.h"
#include"Cell_list.h"
#include"Neigh_list.h"
#include"Initialize_atom_velocity.h"
#include"rescale_atom_velocity.h"
#include"force_calculation.h"
#include"velocity_verlet.h"
#include"Euler_forward.h"
#include"Leap_frog.h"
#include"temperature.h"
#include"strain.h"


int main()
{
//int nx,ny,nz;
//double a;
FILE *read;    //File to read inputs
//Initialisation of variables required
int nx,ny,nz,Crystal_type,NDX,NDY,NDZ,Totaldomains,Totalatoms,t,BC,incremental_time;
double alattice,c,x,y,z,dx,dy,dz,TDesired,Dt,D,alpha,mass,Kb,strainn;
int i,Totaltime,NSave,Integrator;
read = fopen ( "Input.txt","r");

fscanf(read,"%d %d %d %d %lf %lf %d %d %d %d %lf %d %lf %lf %d %d  ",&Crystal_type,&nx,&ny,&nz,&alattice,&c,&NDX,&NDY,&NDZ,&BC,&TDesired,&Totaltime,&D,&alpha,&NSave,&Integrator);
Dt=0.001;
Kb=8.617E-5;
mass=27;   //Aluminium
//printf("check %lg %lf %lf \n",Dt,Kb,mass);   
//BCC
if(Crystal_type==1)             // total number of atoms 
{
Totalatoms=(2*nx*ny*nz);
t=2;
}
//FCC and HCP
else
{
Totalatoms=(4*nx*ny*nz);
t=4;
}
if((Crystal_type==1)||(Crystal_type==2))
{
x=(nx*alattice);				// length in the x direction
y=(ny*alattice);				// length in the Y direction
z=(nz*alattice);				// length in the Z direction

dx=(x)/NDX;				// length of each sub-domain in the x direction
dy=(y)/NDY;				// length of each sub-domain in the Y direction
dz=(z)/NDZ;				// length of each sub-domain in the Z direction

Totaldomains= (NDX*NDY*NDZ);		// total number of sub-domains 

}
else
{
x=(2*nx*alattice);				// length in the x direction
y=((ny)*(sqrt(3))*(alattice));				// length in the Y direction
z=(nz*c);				// length in the Z direction

dx=(x)/NDX;				// length of each sub-domain in the x direction
dy=(y)/NDY;				// length of each sub-domain in the Y direction
dz=(z)/NDZ;				// length of each sub-domain in the Z direction

Totaldomains= (NDX*NDY*NDZ);		// total number of sub-domains 

}

double P[t][3],X[Totalatoms],Y[Totalatoms],Z[Totalatoms];


double Cell_Atoms[Totaldomains][Totalatoms];
int Cell_Num[Totaldomains],Atom_list[Totaldomains][Totalatoms],Atoms_Cell[Totalatoms],Neigh[Totaldomains][27];

double V_rand[Totalatoms],Vx_rand[Totalatoms],Vy_rand[Totalatoms],Vz_rand[Totalatoms],V_scaled,instantaneous_temp,Vx_Scaled[Totalatoms],Vy_Scaled[Totalatoms],Vz_Scaled[Totalatoms],V_sum,Scale_factor;

double Fx[Totalatoms],Fy[Totalatoms],Fz[Totalatoms],Fx_ij[Totalatoms],Fy_ij[Totalatoms],Fz_ij[Totalatoms];
double R_ij,X_ij[Totalatoms][Totalatoms],Y_ij[Totalatoms][Totalatoms],Z_ij[Totalatoms][Totalatoms];
double X_New[Totalatoms],Y_New[Totalatoms],Z_New[Totalatoms],X_Force[Totalatoms],Y_Force[Totalatoms],Z_Force[Totalatoms];
double Velocity,kinetic_energy,Vx_New[Totalatoms],Vy_New[Totalatoms],Vz_New[Totalatoms],X_Old[Totalatoms],Y_Old[Totalatoms],Z_Old[Totalatoms],Potential_energy,PE,Total;
double Old_Vel_X[Totalatoms],Old_Vel_Y[Totalatoms],Old_Vel_Z[Totalatoms],New_Vel_X[Totalatoms],New_Vel_Y[Totalatoms],New_Vel_Z[Totalatoms],Force_y;
//for box
double b1,b2,b3,d1,d2,d3,a1,a2,stress,s,Vy,temp;
double boxs[Totalatoms];
              
        
        s=((alattice*nx)*(alattice*nz));  //surface area

                         
Force_y=0;

Velocity=0;

for(i=0;i<Totalatoms;i++)
{

X_Force[i]=0;
Y_Force[i]=0;
Z_Force[i]=0;


Vx_New[i]=0.0;
Vy_New[i]=0.0;
Vz_New[i]=0.0;
}


//Output File
FILE *output;
output=fopen("velocity_verlet.txt","w");

FILE *output1;
output1=fopen("Total_energy.txt","w");

FILE *output2;
output2=fopen("Euler_forward.txt","w");

FILE *output3;
output3=fopen("Verlet.txt","w");

FILE *output4;
output4=fopen("Leap_frog.txt","w");

FILE *output5;
output5=fopen("Stress_time.txt","w");

FILE *output6;
output6=fopen("Temperature_time.txt","w");

FILE *output7;
output7=fopen("KE.txt","w");

FILE *output8;
output8=fopen("PE.txt","w");

if(Crystal_type==1)
  {
bcc_structure(nx,ny,nz,alattice,Totalatoms,t,P,X,Y,Z);
//Atoms_cell=(Totalatoms/Totaldomains);

//Cell_list(nx,ny,nz,dx,dy,dz,Totaldomains,NDX,NDY,NDZ,Totalatoms,Cell_Atoms,X,Y,Z,Cell_Num,Atom_list,Atoms_Cell,alattice,x,y,z);
//Neigh_list(nx,ny,nz,NDX,NDY,NDZ,Totaldomains,Totalatoms,Neigh,BC);

//New addition
Initialize_atom_velocity(Totalatoms,V_rand,Vx_rand,Vy_rand,Vz_rand);
rescale_atom_velocity(Totalatoms,Vx_rand,Vy_rand,Vz_rand,Kb,mass,TDesired,Vx_Scaled,Vy_Scaled,Vz_Scaled,&V_scaled,V_sum,Scale_factor);
//instantaneous_temperature(Totalatoms,mass,V_scaled,Kb,instantaneous_temp);

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);


 }  
 else if(Crystal_type==2)
    {
fcc_structure(alattice,nx,ny,nz,Totalatoms,t,P,X,Y,Z,boxs,&d1,&d2,&d3);

//Atoms_cell=(Totalatoms/Totaldomains);

//Cell_list(nx,ny,nz,dx,dy,dz,Totaldomains,NDX,NDY,NDZ,Totalatoms,Cell_Atoms,X,Y,Z,Cell_Num,Atom_list,Atoms_Cell,alattice,x,y,z);
//Neigh_list(nx,ny,nz,NDX,NDY,NDZ,Totaldomains,Totalatoms,Neigh,BC);

//New addition
Initialize_atom_velocity(Totalatoms,V_rand,Vx_rand,Vy_rand,Vz_rand);
rescale_atom_velocity(Totalatoms,Vx_rand,Vy_rand,Vz_rand,Kb,mass,TDesired,Vx_Scaled,Vy_Scaled,Vz_Scaled,&V_scaled,V_sum,Scale_factor);
//instantaneous_temperature(Totalatoms,mass,V_scaled,Kb,instantaneous_temp);

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);


    }
 else if(Crystal_type==3)
   { 
hcp_structure(alattice,c,nx,ny,nz,Totalatoms,t,P,X,Y,Z);

//Atoms_cell=(Totalatoms/Totaldomains);
//Cell_list(nx,ny,nz,dx,dy,dz,Totaldomains,NDX,NDY,NDZ,Totalatoms,Cell_Atoms,X,Y,Z,Cell_Num,Atom_list,Atoms_Cell,alattice,x,y,z);
//Neigh_list(nx,ny,nz,NDX,NDY,NDZ,Totaldomains,Totalatoms,Neigh,BC);

//New addition
Initialize_atom_velocity(Totalatoms,V_rand,Vx_rand,Vy_rand,Vz_rand);
rescale_atom_velocity(Totalatoms,Vx_rand,Vy_rand,Vz_rand,Kb,mass,TDesired,Vx_Scaled,Vy_Scaled,Vz_Scaled,&V_scaled,V_sum,Scale_factor);
//instantaneous_temperature(Totalatoms,mass,V_scaled,Kb,instantaneous_temp);

force_calculation(alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy);

   }

for(i=0;i<Totalatoms;i++)
{
printf("Boxs= %lf\n",boxs[i]);

}
printf("d1= %lf d2=%lf d3=%lf \n",d1,d2,d3);

// for creating output file 
int T;
for(T=1;T<=Totaltime;T++)
{
if(Integrator==1)
{
velocity_verlet(alattice,nx,ny,nz,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,X_New,Y_New,Z_New,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,incremental_time,X_Force,Y_Force,Z_Force,Totalatoms,Vx_Scaled,Vy_Scaled,Vz_Scaled,Dt,mass,R_ij,X_ij,Y_ij,Z_ij,output,NSave,T,Vx_New,Vy_New,Vz_New,V_scaled,d1,d2,d3,boxs,Potential_energy);

// to calculate keinetic energy
	for(i=0;i<Totalatoms;i++)
	{
	Velocity=Velocity+(Vx_Scaled[i]*Vx_Scaled[i])+(Vy_Scaled[i]*Vy_Scaled[i])+(Vz_Scaled[i]*Vz_Scaled[i]);
	}
//printf("Velocity %lf\n",Velocity);

kinetic_energy=0.5*mass*Velocity;

fprintf(output7,"\n %d\t %lf\n",T,kinetic_energy);

for(i=0;i<Totalatoms;i++)
{
X[i]=X_New[i];
Y[i]=Y_New[i];
Z[i]=Z_New[i];

}

Potential(R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,output8,T);
//PE=Potential_energy;
//fprintf(output8,"\n %d\t %lf\n",t,Potential_energy);

//printf("Potential Energy : %lf\n",PE);
Total_energy(Total,kinetic_energy,R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,output1,T,d1,d2,d3,boxs,output8);

//Stress calculation
stress_calculation(d1,d2,d3,Fx_ij,Fy_ij,Fz_ij,Force_y,stress,output5,T,Totalatoms,boxs,s);
strain(Totalatoms,alattice,nz,Vy_Scaled,NSave,Dt,Vy,strainn);
temperature(Totalatoms,mass,V_scaled,Kb,temp,T,output6);

}

else if(Integrator==2)
{
Euler_forward(output2,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,incremental_time,mass,Dt,D,alpha,X,Y,Z,X_New,Y_New,Z_New,Vx_Scaled,Vy_Scaled,Vz_Scaled,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,R_ij,X_ij,Y_ij,Z_ij,Vx_New,Vy_New,Vz_New,T,d1,d2,d3,boxs,Potential_energy);

// to calculate keinetic energy
	for(i=0;i<Totalatoms;i++)
	{
	Velocity=Velocity+(Vx_Scaled[i]*Vx_Scaled[i])+(Vy_Scaled[i]*Vy_Scaled[i])+(Vz_Scaled[i]*Vz_Scaled[i]);
        
	}
//printf("Velocity %lf\n",Velocity);

kinetic_energy=0.5*mass*Velocity;

fprintf(output7,"\n %d\t %lf\n",t,kinetic_energy);

for(i=0;i<Totalatoms;i++)
{
X[i]=X_New[i];
Y[i]=Y_New[i];
Z[i]=Z_New[i];

}

Potential(R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,output8,T);
//PE=Potential_energy;
fprintf(output8,"\n %d\t %lf\n",T,Potential_energy);

//printf("Potential Energy : %lf\n",PE);
Total_energy(Total,kinetic_energy,R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,output1,T,d1,d2,d3,boxs,output8);

//Stress calculation
stress_calculation(d1,d2,d3,Fx_ij,Fy_ij,Fz_ij,Force_y,stress,output5,T,Totalatoms,boxs,s);
strain(Totalatoms,alattice,nz,Vy_Scaled,NSave,Dt,Vy,strainn);
temperature(Totalatoms,mass,V_scaled,Kb,temp,T,output6);

}
else if(Integrator==3)
{
Verlet(alattice,nx,ny,nz,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,X_New,Y_New,Z_New,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,incremental_time,X_Force,Y_Force,Z_Force,Totalatoms,Vx_Scaled,Vy_Scaled,Vz_Scaled,Dt,mass,R_ij,X_ij,Y_ij,Z_ij,output,NSave,T,Vx_New,Vy_New,Vz_New,X_Old,Y_Old,Z_Old,output3,d1,d2,d3,boxs,Potential_energy);

// to calculate keinetic energy
	for(i=0;i<Totalatoms;i++)
	{
	Velocity=Velocity+(Vx_Scaled[i]*Vx_Scaled[i])+(Vy_Scaled[i]*Vy_Scaled[i])+(Vz_Scaled[i]*Vz_Scaled[i]);
	}
//printf("Velocity %lf\n",Velocity);

kinetic_energy=0.5*mass*Velocity;

fprintf(output7,"\n %d\t %lf\n",t,kinetic_energy);

for(i=0;i<Totalatoms;i++)
{
X[i]=X_New[i];
Y[i]=Y_New[i];
Z[i]=Z_New[i];

}

Potential(R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,output8,T);
//PE=Potential_energy;
//fprintf(output8,"\n %d\t %lf\n",t,Potential_energy);

//printf("Potential Energy : %lf\n",PE);
Total_energy(Total,kinetic_energy,R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,output1,T,d1,d2,d3,boxs,output8);

//Stress calculation
stress_calculation(d1,d2,d3,Fx_ij,Fy_ij,Fz_ij,Force_y,stress,output5,T,Totalatoms,boxs,s);
strain(Totalatoms,alattice,nz,Vy_Scaled,NSave,Dt,Vy,strainn);
temperature(Totalatoms,mass,V_scaled,Kb,temp,T,output6);

}

else if(Integrator==4)
{
Leap_frog(alattice,nx,ny,nz,NDX,NDY,NDZ,Totaldomains,D,alpha,X,Y,Z,X_New,Y_New,Z_New,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,incremental_time,X_Force,Y_Force,Z_Force,Totalatoms,Vx_Scaled,Vy_Scaled,Vz_Scaled,Dt,mass,R_ij,X_ij,Y_ij,Z_ij,NSave,T,Vx_New,Vy_New,Vz_New,X_Old,Y_Old,Z_Old,output4,Old_Vel_X,Old_Vel_Y,Old_Vel_Z,New_Vel_X,New_Vel_Y,New_Vel_Z,d1,d2,d3,boxs,Potential_energy);

// to calculate keinetic energy
	for(i=0;i<Totalatoms;i++)
	{
	Velocity=Velocity+(Vx_Scaled[i]*Vx_Scaled[i])+(Vy_Scaled[i]*Vy_Scaled[i])+(Vz_Scaled[i]*Vz_Scaled[i]);
	}
//printf("Velocity %lf\n",Velocity);

kinetic_energy=0.5*mass*Velocity;

fprintf(output7,"\n %d\t %lf\n",t,kinetic_energy);

for(i=0;i<Totalatoms;i++)
{
X[i]=X_New[i];
Y[i]=Y_New[i];
Z[i]=Z_New[i];

}

Potential(R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,d1,d2,d3,boxs,Potential_energy,T);
//PE=Potential_energy;

//fprintf(output8,"\n %d\t %lf\n",t,Potential_energy);

//printf("Potential Energy : %lf\n",PE);
Total_energy(Total,kinetic_energy,R_ij,alpha,D,Potential_energy,alattice,nx,ny,nz,Totalatoms,NDX,NDY,NDZ,Totaldomains,X,Y,Z,Fx,Fy,Fz,Fx_ij,Fy_ij,Fz_ij,X_ij,Y_ij,Z_ij,output1,T,d1,d2,d3,boxs,output8);

//Stress calculation
stress_calculation(d1,d2,d3,Fx_ij,Fy_ij,Fz_ij,Force_y,stress,output5,T,Totalatoms,boxs,s);
strain(Totalatoms,alattice,nz,Vy_Scaled,NSave,Dt,Vy,strainn);
temperature(Totalatoms,mass,V_scaled,Kb,temp,T,output6);

}




}


return 0;
}
