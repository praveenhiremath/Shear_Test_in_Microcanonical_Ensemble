#ifndef TOTAL_ENERGY_H_INCLUDED
#define TOTAL_ENERGY_H_INCLUDED

void Total_energy(double Total,double kinetic_energy,double R_ij,double alpha,double D,double Potential_energy,double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,double X[],double Y[],double Z[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],,FILE *output1,int T,double d1,double d2,double d3,double boxs[],FILE *output8);

#endif
