
#ifndef REPULSIVE_BROWNIAN_H
#define REPULSIVE_BROWNIAN_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"
#include "brownian.h"

#include <omp.h>

namespace REPULSIVE_BROWNIAN
{
  struct
  TEMPORAL_VARIABLE; // tell compiler to expect struct def.
  
  long
  main_EQUILIBRATION(TRAJECTORY& TRAJ,
		     POTENTIAL_SET& POTs,
		     RECORD_DATA& DATA,
		     COND& given_condition);

  double
  record_simulation_data(RECORD_DATA& DATA,
			 TRAJECTORY& TRAJ, const long index_t_now,
			 MATRIX& energy);

  double
  report_simulation_info(TRAJECTORY& TRAJ,
			 MATRIX& energy,
			 TEMPORAL_VARIABLE& VAR);
  double
  OMP_compute_RDIST(TRAJECTORY& TRAJ, const long index_t_now,
		    RDIST& R_boost, long* tmp_index_vec,
		    const long N_THREADS_BD);
  double
  OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const long index_t_now, const long index_t_next,
			   POTENTIAL_SET& POTs, MATRIX* force_repulsion, MATRIX* force_random,
			   RDIST& R_boost, MATRIX* vec_boost_Nd_parallel,
			   RNG_BOOST& RNG,
			   const long N_THREADS_BD,
			   COND& given_condition, TEMPORAL_VARIABLE& VAR);

  // inline functions
  double
  sum_virial_components(MATRIX& energy);
}

struct REPULSIVE_BROWNIAN::
  TEMPORAL_VARIABLE
  : BROWNIAN::BROWNIAN_VARIABLE
{
  // the variables on here is supporting for interface
  // do not use the core variable inside this structure
  // for designing purpose, it is of importance to use individual variable definition outside of this structure
  MATRIX *vec_boost_Nd_parallel;
  MATRIX *force_repulsion;

  long N_skip_rdist = 0;
  
  double RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz;
  double RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz;

  double energy_repulsive_potential;
  
  // member functions for virials
  long
  virial_initial();
  double
  record_virial_into_energy_array
  (MATRIX& energy);
  
  // related with shear flow
  double Wi_tau_R;
  TEMPORAL_VARIABLE()
  {
    std::cout << "ERR: Basic constructor for TEMPORAL_VARIABLE structure is not supported\n";
  }
  TEMPORAL_VARIABLE
  (COND& given_condition, long given_N_basic);
  virtual
  ~TEMPORAL_VARIABLE();

};

inline long
REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
virial_initial()
{
  // note that this function is not inherite from previous function, directly
  // rather than, it explicitly call the function defined in the function in parents class
  BROWNIAN::BROWNIAN_VARIABLE::
    virial_initial();
  
  RF_repulsion_xx = 0.; RF_repulsion_yy = 0.; RF_repulsion_zz = 0.;
  RF_repulsion_xy = 0.; RF_repulsion_xz = 0.; RF_repulsion_yz = 0.;

  energy_repulsive_potential = 0.;
  
  return 0;
}

inline double
REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
record_virial_into_energy_array(MATRIX& energy)
{
  double time_st = omp_get_wtime();

  BROWNIAN::BROWNIAN_VARIABLE::
    record_virial_into_energy_array(energy);
  
  energy(18) = RF_repulsion_xx;
  energy(19) = RF_repulsion_yy;
  energy(20) = RF_repulsion_zz;

  energy(21) = RF_repulsion_xy;
  energy(22) = RF_repulsion_xz;
  energy(23) = RF_repulsion_yz;

  energy(2) = energy_repulsive_potential;
  
  return omp_get_wtime() - time_st;
}
				
inline double
REPULSIVE_BROWNIAN::
sum_virial_components(MATRIX& energy)
{
  double time_st = omp_get_wtime();
  BROWNIAN::sum_virial_components(energy);
  
  for(long i=6; i<12; i++)
    {
      energy(i) += energy(i + 6*2);
    }


  energy(1) = energy(2) + energy(3);
  
  return omp_get_wtime() - time_st;
}

#endif
