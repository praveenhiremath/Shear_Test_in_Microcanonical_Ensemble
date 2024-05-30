#ifndef CELL_LIST_H_INCLUDED
#define CELL_LIST_H_INCLUDED
void Cell_list(int nx,int ny,int nz,double dx,double dy,double dz,int Totaldomains,int NDX,int NDY,int NDZ,int Totalatoms,double Cell_Atoms[Totaldomains][Totalatoms],double X[],double Y[],double Z[],int Cell_Num[],int Atom_list[Totaldomains][Totalatoms],int Atoms_Cell[],double alattice,double x,double y,double z);


#endif
