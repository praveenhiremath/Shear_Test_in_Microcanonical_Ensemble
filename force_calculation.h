#ifndef FORCE_CALCULATION_H_INCLUDED
#define FORCE_CALCULATION_H_INCLUDED


void force_calculation(double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],double d1,double d2,double d3,double boxs[],double Potential_energy);

#endif
