#include<stdio.h>
#include<math.h>
#include<stdlib.h>
//#include"rescale_atom_velocity.h"

void rescale_atom_velocity(int Totalatoms, double Vx_rand[],double Vy_rand[],double Vz_rand[],double Kb, double mass, int TDesired,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double *V_scaled,double V_sum,double Scale_factor)

{
int i;
double Vx_sum,Vy_sum,Vz_sum,Vx_sum_atom,Vy_sum_atom,Vz_sum_atom;
double lmx,lmy,lmz,instantaneous_temp;
lmx=0.0; lmy=0.0; lmz=0.0;

Vx_sum =0;
Vy_sum =0;
Vz_sum =0;
V_sum  =0;

for(i=0;i<Totalatoms;i++)
{
//Sum of all random velocities
Vx_sum=Vx_sum+Vx_rand[i];
Vy_sum=Vy_sum+Vy_rand[i];
Vz_sum=Vz_sum+Vz_rand[i];

//Sum of all squared random velocities
V_sum =V_sum+(Vx_rand[i]*Vx_rand[i])+(Vy_rand[i]*Vy_rand[i])+(Vz_rand[i]*Vz_rand[i]);

}

//Mean velocity of an atom in X,Y,Z directions		
	Vx_sum_atom = Vx_sum/Totalatoms;
	Vy_sum_atom = Vy_sum/Totalatoms;
	Vz_sum_atom = Vz_sum/Totalatoms;

//Scaling factor
	Scale_factor = sqrt((3*Kb*Totalatoms*TDesired)/(V_sum*mass));
        //Scale_factor = sqrt((3*TDesired)/(V_sum));
        printf("Scaling factor is %lf\n",Scale_factor);

//Scaled velocity componenets
	for(i=0;i<Totalatoms;i++)
	{				
	   Vx_Scaled[i]=(Vx_rand[i]-Vx_sum_atom)*Scale_factor; 
	   Vy_Scaled[i]=(Vy_rand[i]-Vy_sum_atom)*Scale_factor;
	   Vz_Scaled[i]=(Vz_rand[i]-Vz_sum_atom)*Scale_factor;
          //printf("Scaled Velocity in X %lf \n",Vx_Scaled[i]);
          //printf("Scaled Velocity in Y %lf \n",Vy_Scaled[i]);
          //printf("Scaled Velocity in Z %lf \n",Vz_Scaled[i]);

     lmx=lmx+Vx_Scaled[i];
     lmy=lmy+Vy_Scaled[i];
     lmz=lmz+Vz_Scaled[i];

	}



printf(" Linear momentums: %lf %lf %lf\n",lmx,lmy,lmz);

*V_scaled=0.0;
for(i = 0 ; i <Totalatoms ; i++)
		    {
		        //printf("%d %lf %lf %lf\n",i+1, Vx_Scaled[i], Vy_Scaled[i], Vz_Scaled[i]);
//Sum of sqaured velocities of atoms                         
                        *V_scaled=*V_scaled+(Vx_Scaled[i]*Vx_Scaled[i]+Vy_Scaled[i]*Vy_Scaled[i]+Vz_Scaled[i]*Vz_Scaled[i]);
//                      printf("Scaled Velocity is %lf\n",V_scaled);
		    }

//printf("Scaled Velocity is %lf\n",V_scaled);


}


