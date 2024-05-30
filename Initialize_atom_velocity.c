#include<stdio.h>
#include<math.h>
#include<stdlib.h>


void Initialize_atom_velocity(int Totalatoms,double V_rand[],double Vx_rand[],double Vy_rand[],double Vz_rand[])
{
int i;
//double V_rand[],Vx_rand[],Vy_rand[],Vz_rand[];
double max,min;
max=0.5; 
min=-0.5;

for(i=0;i<Totalatoms;i++)
{
Vx_rand[i] =rand()*(max-min)/RAND_MAX+min;
Vy_rand[i] =rand()*(max-min)/RAND_MAX+min;
Vz_rand[i] =rand()*(max-min)/RAND_MAX+min;


}

V_rand[Totalatoms]=0.0;
for(i=0;i<Totalatoms;i++)
{
V_rand[i]=V_rand[i]+(Vx_rand[i]*Vx_rand[i]+Vy_rand[i]*Vy_rand[i]+Vz_rand[i]*Vz_rand[i]);
//printf("Random velocity  %lf\n",V_rand[i]);

}

}
