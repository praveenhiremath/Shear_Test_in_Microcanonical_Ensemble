#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include"Cell_list.c"
//#include"Neigh_list.c"


void force_calculation(double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],double d1,double d2,double d3,double boxs[],double Potential_energy)
{

//double R_ij[Totalatoms][Totalatoms],X_ij[Totalatoms][Totalatoms],Y_ij[Totalatoms][Totalatoms],Z_ij[Totalatoms][Totalatoms];
int i,k,l,j,s;
double R_0,C1,C2,lx,ly,lz;
double A1,A2;
R_0=3.4068;


//Force on 'i'th atom due to 'j'th atoms
//Fx_ij[0]=0.0;
//Fy_ij[0]=0.0;
//Fz_ij[0]=0.0;

lx=alattice*nx;
ly=alattice*ny;
lz=alattice*nz;
//printf("lx,ly and lz %lf %lf %lf\n", lx,ly,lz);

//To find distance between any two atoms in terms of "r" and also in terms of their components
//force on 'i'th atom due to 'j'th atom 


for(i=0;i<Totalatoms-1;i++)
{
  for(j=i+1;j<Totalatoms;j++)
  {
    X_ij[i][j]=X[i]-X[j];
    Y_ij[i][j]=Y[i]-Y[j];
    Z_ij[i][j]=Z[i]-Z[j];
//    R_ij[i][j]=sqrt((X_ij[i][j]*X_ij[i][j])+(Y_ij[i][j]*Y_ij[i][j])+(Z_ij[i][j]*Z_ij[i][j]));

   if(fabs(X_ij[i][j])>(lx/2))
    { 
      if(X_ij[i][j]>0)
//      {
      lx=lx;
//      }
      else if(X_ij[i][j]<0)
//      {
       lx=(-lx);
//       }
      
      X_ij[i][j]=X_ij[i][j]-lx;
    }
   if(fabs(Y_ij[i][j])>(ly/2))
    { 
      if(Y_ij[i][j]>0)
//      {
      ly=ly;
//      }
      else if(Y_ij[i][j]<0)
//      {
       ly=(-ly);
//       }
      
      Y_ij[i][j]=Y_ij[i][j]-ly;
    }
   
    R_ij=sqrt((X_ij[i][j]*X_ij[i][j])+(Y_ij[i][j]*Y_ij[i][j])+(Z_ij[i][j]*Z_ij[i][j]));
//printf("R_ij : %lf\n",R_ij[i][j]);
//printf("Atomic distances : %lf %lf %lf\n",X_ij[i][j],Y_ij[i][j],Z_ij[i][j]);

//Morse' force calculation
//     if(R_ij<R_0)
//     {
   C1=exp(-2*alpha*(R_ij-R_0));
   C2=exp(-alpha*(R_ij-R_0));
     Fx_ij[i]=Fx_ij[i]+(D*((-2*alpha*C1)+(2*alpha*C2)))*(X_ij[i][j]/R_ij);
     Fy_ij[i]=Fy_ij[i]+(D*((-2*alpha*C1)+(2*alpha*C2)))*(Y_ij[i][j]/R_ij);
     Fz_ij[i]=Fz_ij[i]+(D*((-2*alpha*C1)+(2*alpha*C2)))*(Z_ij[i][j]/R_ij);


//Potential energy
A1=exp((-2*alpha)*(R_ij-R_0));
A2=exp(((-alpha)*(R_ij-R_0)));
Potential_energy= (D)*(A1-(2*A2));
//printf("Potential is %lf\n",Potential_energy);   //Potential of each atom

//printf("Morse forces on %d th atom: %lf %lf %lf \n",i,Fx_ij[i],Fy_ij[i],Fz_ij[i]);


     Fx_ij[j]=Fx_ij[j]-(D*((-2*alpha*C1)+(2*alpha*C2)))*(X_ij[j][j+1]/R_ij);
     Fy_ij[j]=Fy_ij[j]-(D*((-2*alpha*C1)+(2*alpha*C2)))*(Y_ij[j][j+1]/R_ij);
     Fz_ij[j]=Fz_ij[j]-(D*((-2*alpha*C1)+(2*alpha*C2)))*(Z_ij[j][j+1]/R_ij);

//printf("Morse forces on %d th atom: %lf %lf %lf \n",j,Fx_ij[j],Fy_ij[j],Fz_ij[j]);

//   }

  }
}


//Bottom region

for(s=0;s<d1;s++)
{
i=boxs[s];
  
     Fx_ij[i]=0.0;
     Fy_ij[i]=0.0;
     Fz_ij[i]=0.0;

 }

//Middle region

for(s=d1;s<d2;s++)
{
i=boxs[s];
  
     Fx_ij[i]=0.0;
     Fy_ij[i]=0.0;
     Fz_ij[i]=0.0;

}


//Top region
for(s=d2;s<d3;s++)
{
i=boxs[s];
 
     Fx_ij[i]=0.0;
     Fy_ij[i]=Fy_ij[i];
     Fz_ij[i]=0.0;


}

}
