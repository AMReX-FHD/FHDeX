#ifndef INS_functions_F_H_
#define _INS_functions_F_H_

#ifdef __cplusplus

extern "C" 
{
#endif

// in membrane_forces.F90


   void user_force_calc_init(const char* inputs_file, const int length);

   void user_force_calc(const amrex_real* spec3xPos, const amrex_real* spec3yPos, const amrex_real* spec3zPos, amrex_real* spec3xForce, amrex_real* spec3yForce, amrex_real* spec3zForce, const int* length, const int* step, const species* particleInfo);

   void user_force_calc_destroy();


#ifdef __cplusplus
}
#endif

#endif
