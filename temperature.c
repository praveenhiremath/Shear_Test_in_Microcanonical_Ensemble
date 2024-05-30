#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//temperature 
void temperature(int Totalatoms,double mass,double V_scaled,double Kb,double temp,int T,FILE *output6)
{

temp=(mass*V_scaled)/(3*Totalatoms*Kb);
fprintf(output6," %d %lf \n",T,temp);
//printf(" %d %lf \n",T,temp);
}

