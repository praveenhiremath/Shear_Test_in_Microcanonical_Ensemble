#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void stress_calculation(double d1,double d2,double d3,double Fx_ij[],double Fy_ij[],double Fz_ij[],double Force_y,double stress,FILE *output5,int T,int Totalatoms,double boxs[],int s)
{
int i,p;


for(p=d2+1;p<=d3;p++)
{
i=boxs[p];
printf("force1 %lf\n",Fy_ij[i]);
//for(i=d2;i<d3;i++)


Force_y=Fy_ij[i]+Force_y;
//printf("force %lf\n",Force_y);

}
//stress calculation
stress=(Force_y/s);
//printf("Stress: %lf\n",stress);
fprintf(output5," %d   %lf\n",T,stress);

}
