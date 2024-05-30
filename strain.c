#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"strain.h"

void strain(int Totalatoms,double alattice,int nz,double Vy_Scaled[],double NSave,double Dt,double Vy,double strainn)
{
double lz;
int i;
lz=nz*alattice;
Vy=0;
for(i=0;i<Totalatoms;i++)
{
Vy=Vy+Vy_Scaled[i];

}

strainn=(Vy*NSave*Dt)/lz;
printf("strain %lf  \n",strainn);

}
