//BCC crystal structure

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void bcc_structure(int nx,int ny,int nz,double alattice,int Totalatoms,int t,double P[t][3],double X[],double Y[],double Z[])

{
FILE *bcc;
bcc = fopen ( "bcc.txt", "w");   //output file containing BCC crystal structure,without cell data
int i,j,k,l,R;

//Ovito format
fprintf(bcc, "Position data for file\n\n%d atoms\n",2*nx*ny*nz);
fprintf(bcc,"\n1 atom types\n");
fprintf(bcc,"0.000000   %lf xlo xhi\n", alattice*nx);
fprintf(bcc,"0.000000   %lf ylo yhi\n", alattice*ny);
fprintf(bcc,"0.000000   %lf zlo zhi\n", alattice*nz);
fprintf(bcc,"\nAtoms\n\n");


//motif atom co-ordinates
P[0][0]=0,P[0][1]=0,P[0][2]=0,P[1][0]=alattice/2,P[1][1]=alattice/2,P[1][2]=alattice/2;
R=0;
l=1;
for (i=0;i<nz;i++)
{
    for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<2;t++)
         {
          
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*alattice);
          Z[R]=P[t][2]+(i*alattice);
fprintf(bcc,"%d 1 %f %f %f\n",l,X[R],Y[R],Z[R]);
l++;         
R++;
          }
        }
     }
}
 
/*printf("Totalatoms: %d \n",Totalatoms);
for(i=0;i<Totalatoms;i++)
 {                              printf("Px[i]:%d %lf \n",i,X[i]);
                               printf("Py[y]:%d %lf \n",i,Y[i]);
                               printf("Pz[z]:%d %lf \n",i,Z[i]);
 }
 */      
}
