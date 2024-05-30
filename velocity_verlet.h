#ifndef VELOCITY_VERLET_H_INCLUDED
#define VELOCITY_VERLET_H_INCLUDED


void velocity_verlet(double alattice,int nx, int ny, int nz,int NDX,int NDY,int NDZ,int Totaldomains,double D,double alpha,double X[],double Y[],double Z[],double X_New[],double Y_New[],double Z_New[],double Fx[],double Fy[],double Fz[],double Fx_ij[],double Fy_ij[],double Fz_ij[],int incremental_time,double X_Force[],double Y_Force[],double Z_Force[],int Totalatoms,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double Dt,double mass,double R_ij,double X_ij[Totalatoms][Totalatoms],double Y_ij[Totalatoms][Totalatoms],double Z_ij[Totalatoms][Totalatoms],FILE *output,int NSave,int T,double Vx_New[],double Vy_New[],double Vz_New[],double V_scaled,double d1,double d2,double d3,double boxs[],double Potential_energy);

#endif
