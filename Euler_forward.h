#ifndef EULER_FORWARD_H_INCLUDED
#define EULER_FORWARD_H_INCLUDED

void Euler_forward(FILE *output2,double alattice,int nx, int ny, int nz,int Totalatoms,int NDX,int NDY,int NDZ,int Totaldomains,int incremental_time,double mass,double Dt,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],double Vx_New[],double Vy_New[],double Vz_New[],int T,double d1,double d2,double d3,double boxs[],double Potential_energy);

#endif
