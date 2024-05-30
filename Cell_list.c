#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void Cell_list(int nx,int ny,int nz,double dx,double dy,double dz,int Totaldomains,int NDX,int NDY,int NDZ,int Totalatoms,double Cell_Atoms[Totaldomains][Totalatoms],double X[],double Y[],double Z[],int Cell_Num[],int Atom_list[Totaldomains][Totalatoms],int Atoms_Cell[],double alattice,double x,double y,double z)
{
//double Cell_Neigh[Totaldomains][27],global_id;
//double Cell[NDX][NDY][NDZ];

int i,j,k,E,Ix,Iy,Iz;
FILE *out;
out = fopen ( "subdomain.txt", "w");

//printf("Total atoms %d Total domains %d \n",Totalatoms,Totaldomains);
//OVITO format output
fprintf(out,"ITEM: TIMESTEP\n");
fprintf(out,"0\n");
fprintf(out,"ITEM: NUMBER OF ATOMS\n");
fprintf(out,"%d \n",Totalatoms);
fprintf(out,"ITEM: BOX BOUNDS\n");
fprintf(out,"0.000000   %lf xlo xhi\n", x);
fprintf(out,"0.000000   %lf ylo yhi\n", y);
fprintf(out,"0.000000   %lf zlo zhi\n", z);
fprintf(out,"ITEM: ATOMS id x y z cell\n");


j=0;


     for(i=0;i<Totalatoms;i++)

       
                    {
            			Ix= (int)(X[i]/dx)+1;
       				Iy= (int)(Y[i]/dy)+1;
       				Iz= (int)(Z[i]/dz)+1;
       				E=  Iy + NDY*(Ix-1+NDX*(Iz-1));   //cell_id to which 'i'th atom belongs
                              Cell_Num[i]=E;                      //array containing cell_data
                              //printf("Cell_Num[%d]=%d\n",i,Cell_Num[i]);

  }

//OVITO format output

for(i=0;i<Totalatoms;i++)
{

fprintf(out,"%d %f %f %f %d\n",i+1,X[i],Y[i],Z[i],Cell_Num[i]);

}



//Number of atoms in a cell is given by Atoms_Cell[]

for(j=1; j<=Totaldomains; j++)
		{
		k=0;
			for(i=0; i<Totalatoms; i++)
			{
				if(Cell_Num[i]==j)
				{				
				Atom_list[j][k]=i;     //array representing 'i'th atom belongs to 'j'th cell.
			//printf("In the Cell %d Atom %d is present\n",Cell_Num[i], Atom_list[j][k]);				
			//	printf("%d %d %lf\n",j,k,table1[j][k]);				
				k++;
				}
			}
		
		Atoms_Cell[j]=k;	   //Array representing how many atoms are in 'j'th cell	
		
	//printf("In the cell %d atoms present= %d \n",j,Atoms_Cell[j]);
			}

}
