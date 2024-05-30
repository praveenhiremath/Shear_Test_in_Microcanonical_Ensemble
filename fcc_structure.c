//FCC Crystal structure
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void fcc_structure(double alattice,int nx, int ny, int nz,int Totalatoms,int t,double P[t][3],double X[],double Y[],double Z[],double boxs[],double *d1,double *d2,double *d3)
{
FILE *fcc;
fcc = fopen ( "fcc.txt", "w");   //output file containing FCC crystal structure,without cell data
int i,j,k,l,R,Q;

//Ovito format
fprintf(fcc, "Position data for file\n\n%d atoms\n",4*nx*ny*nz);
fprintf(fcc,"\n1 atom types\n");
fprintf(fcc,"0.000000   %lf xlo xhi\n", alattice*nx);
fprintf(fcc,"0.000000   %lf ylo yhi\n", alattice*ny);
fprintf(fcc,"0.000000   %lf zlo zhi\n", alattice*nz);
fprintf(fcc,"\nAtoms\n\n");



//motif atom co-ordinates
P[0][0]=0,P[0][1]=0,P[0][2]=0,P[1][0]=0,P[1][1]=alattice/2,P[1][2]=alattice/2,P[2][0]=alattice/2,P[2][1]=0,P[2][2]=alattice/2,P[3][0]=alattice/2,P[3][1]=alattice/2,P[3][2]=0;
l=1;
R=0;
for (i=0;i<nz;i++)
{
    for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<4;t++)
           {
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*alattice);
          Z[R]=P[t][2]+(i*alattice);
fprintf(fcc,"%d 1 %f %f %f\n",l,X[R],Y[R],Z[R]);
l++;          
R++;       
           }
         } 
     }
}

Q=0;
//Bottom Region 
for (i=0;i<2;i++)
{
  for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<4;t++)
           {
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*alattice);
          Z[R]=P[t][2]+(i*alattice);

boxs[Q]= R;     
Q++;       
           }
         } 
     }
}
*d1=Q-1;

//Top Region 
for (i=(nz-3);i<nz;i++)
{
  for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<4;t++)
           {
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*alattice);
          Z[R]=P[t][2]+(i*alattice);

boxs[Q]= R;     
Q++;           }
         } 
     }
}
*d2=Q-1;


//Middle Region 
for (i=2;i<(nz-2);i++)
{
  for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<4;t++)
           {
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*alattice);
          Z[R]=P[t][2]+(i*alattice);

boxs[Q]= R;     
Q++;           }
         } 
     }
}
*d3=Q-1;


printf("d1= %lf d2=%lf d3=%lf \n",*d1,*d2,*d3);
}


