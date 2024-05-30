#ifndef RESCALE_ATOM_VELOCITY_H_INCLUDED
#define RESCALE_ATOM_VELOCITY_H_INCLUDED

void rescale_atom_velocity(int Totalatoms, double Vx_rand[],double Vy_rand[],double Vz_rand[],double Kb, double mass, int TDesired,double Vx_Scaled[],double Vy_Scaled[],double Vz_Scaled[],double *V_scaled,double V_sum,double Scale_factor);

#endif
