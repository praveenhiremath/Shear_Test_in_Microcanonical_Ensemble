//HCP crystal structure
#include<stdio.h>
#include<stdlib.h>
#include<math.h>



void hcp_structure(double alattice,double c,int nx,int ny,int nz,int Totalatoms,int t,double P[t][3],double X[],double Y[],double Z[])
{
FILE *hcp;
hcp = fopen ( "hcp.txt", "w");    //output file containing HCP crystal structure,without cell data
int i,j,k,l,R;

//Ovito format
fprintf(hcp, "Position data for file\n\n%d atoms\n",4*nx*ny*nz);
fprintf(hcp,"\n1 atom types\n");
fprintf(hcp,"0.000000   %lf xlo xhi\n", (alattice*nx));
fprintf(hcp,"0.000000   %lf ylo yhi\n", (alattice*(sqrt(3)))*ny);
fprintf(hcp,"0.000000   %lf zlo zhi\n", (c*nz));
fprintf(hcp,"\nAtoms\n\n");


//motif atom co-ordinates
P[0][0]=0,P[0][1]=0,P[0][2]=0,P[1][0]=(alattice/2),P[1][1]=((sqrt(3)*alattice)/2),P[1][2]=0;
P[2][0]=0,P[2][1]=(-alattice/sqrt(3)),P[2][2]=(c/2),P[3][0]=(alattice/2),P[3][1]=(alattice/(2*(sqrt(3)))),P[3][2]=(c/2);

R=0;
l=1;
for (i=0;i<nz;i++)
{
    for  (j=0;j<ny;j++)
    {
        for (k=0;k<nx;k++)
        {
          for(t=0;t<4;t++)
           {
          X[R]=P[t][0]+(k*alattice);
          Y[R]=P[t][1]+(j*((sqrt(3))*alattice));
          Z[R]=P[t][2]+(i*c);
fprintf(hcp,"%d 1 %lf %lf %lf\n",l,X[R],Y[R],Z[R]);
l++;
R++;

           }
        }
     }
}

}
