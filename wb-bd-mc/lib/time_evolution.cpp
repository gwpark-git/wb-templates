
#include "time_evolution.h"

// double INTEGRATOR::time_evolution_Euler(TRAJECTORY& TRAJ, const long index_t_now, const long index_t_next, const long index_bead, MATRIX** contributions, const long N_forces)
// {
//   // contributions[0] represent the Euler
//   // contributions[1] represent mechanical perturbation (if exist)
//   // contributions[2] represent repulsive potential (if exist)
//   // contributions[3] represent interaction through association (if exist)
//   for(long k=0; k<TRAJ.N_dimension; k++)
//     {
//       TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k); // inheritance the current positions
//       for(long ind_f=0; ind_f<N_forces; ind_f++)
//         {
//           // TRAJ(index_t_next, i, k) += TRAJ.dt*force_array[ind_f][i](k);
//         }
//     }
//   return 0;
// }

double
INTEGRATOR::EULER_ASSOCIATION::
cal_connector_force_boost
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      MATRIX& rel_vector = R_minimum_vec_boost[given_index][target_index];
      double distance = R_minimum_distance_boost[given_index](target_index);
      double force = CONNECT.weight[given_index](j)*POTs.f_connector(distance, POTs.force_variables);
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force/distance, // a
                  rel_vector.data, // x
                  1,
                  given_vec.data, // y
                  1);
      
    }
  return omp_get_wtime() - time_st;
}

double
INTEGRATOR::EULER_ASSOCIATION::
cal_connector_force_boost_with_RF
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost,
 double& RF_connector_xx, double& RF_connector_yy, double& RF_connector_zz,
 double& RF_connector_xy, double& RF_connector_xz, double& RF_connector_yz,
 double& energy_elastic_potential)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      MATRIX& rel_vector = R_minimum_vec_boost[given_index][target_index];
      double distance = R_minimum_distance_boost[given_index](target_index);
      double force = CONNECT.weight[given_index](j)*POTs.f_connector(distance, POTs.force_variables);
      double force_div_distance = force/distance;
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force_div_distance, // a
                  rel_vector.data, // x
                  1,
                  given_vec.data, // y
                  1);

      RF_connector_xx += force_div_distance*rel_vector(0)*rel_vector(0);
      RF_connector_yy += force_div_distance*rel_vector(1)*rel_vector(1);
      RF_connector_zz += force_div_distance*rel_vector(2)*rel_vector(2);

      RF_connector_xy += force_div_distance*rel_vector(0)*rel_vector(1);
      RF_connector_xz += force_div_distance*rel_vector(0)*rel_vector(2);
      RF_connector_yz += force_div_distance*rel_vector(1)*rel_vector(2);

      energy_elastic_potential += CONNECT.weight[given_index](j)*POTs.e_connector(distance, POTs.force_variables);
      
    }
  return omp_get_wtime() - time_st;
}



double
INTEGRATOR::EULER::
cal_connector_force_boost
(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& given_vec, long given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      MATRIX& rel_vector = R_minimum_vec_boost[given_index][target_index];
      double distance = R_minimum_distance_boost[given_index](target_index);
      double force = POTs.f_connector(distance, POTs.force_variables);
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force/distance, // a
                  rel_vector.data, // x
                  1,
                  given_vec.data, // y
                  1);
      
    }
  return omp_get_wtime() - time_st;
}


double
INTEGRATOR::EULER::
cal_connector_force_boost_with_RF
(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& given_vec, long given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost,
 double& RF_connector_xx, double& RF_connector_yy, double& RF_connector_zz,
 double& RF_connector_xy, double& RF_connector_xz, double& RF_connector_yz,
 double& energy_elastic_potential)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  for (long j=1; j<CONNECT.TOKEN[given_index]; j++)
    {
      long target_index = (long)CONNECT.HASH[given_index](j);
      MATRIX& rel_vector = R_minimum_vec_boost[given_index][target_index];
      double distance = R_minimum_distance_boost[given_index](target_index);
      double force = POTs.f_connector(distance, POTs.force_variables);
      double force_div_distance = force/distance;
      // a*x + y -> y
      cblas_daxpy(given_vec.size,
                  force_div_distance, // a
                  rel_vector.data, // x
                  1,
                  given_vec.data, // y
                  1);

      RF_connector_xx += force_div_distance*rel_vector(0)*rel_vector(0);
      RF_connector_yy += force_div_distance*rel_vector(1)*rel_vector(1);
      RF_connector_zz += force_div_distance*rel_vector(2)*rel_vector(2);

      RF_connector_xy += force_div_distance*rel_vector(0)*rel_vector(1);
      RF_connector_xz += force_div_distance*rel_vector(0)*rel_vector(2);
      RF_connector_yz += force_div_distance*rel_vector(1)*rel_vector(2);

      energy_elastic_potential += POTs.e_connector(distance, POTs.force_variables);
      
    }
  return omp_get_wtime() - time_st;
}


double
INTEGRATOR::EULER::
cal_repulsion_force_R_boost
(POTENTIAL_SET& POTs, MATRIX& given_vec, long index_particle, RDIST& R_boost)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  long cell_index_particle = R_boost.cell_index[index_particle];
  for(long k=0; k<R_boost.N_neighbor_cells; k++)
    {
      long cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(long p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          long index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.Rsca[index_particle](index_target);
          if (index_target != index_particle)
            {
              // double dimensionless_C_rep =
              double repulsion = POTs.f_repulsion_coupling(POTs.force_variables, index_particle, index_target)*POTs.f_repulsion(distance, POTs.force_variables);
              // double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
	      
              // for(long d=0; d<R_boost.N_dimension; d++)
              //   {
              //     given_vec(d) += repulsion*R_boost.Rvec[index_particle][index_target](d);
              //   }
	      // the following make overhead.
              cblas_daxpy(given_vec.size,
                          repulsion/distance,
                          R_boost.Rvec[index_particle][index_target].data,
                          1,
                          given_vec.data,
                          1);
            }
        }
    }
  return omp_get_wtime() - time_st;
}


double
INTEGRATOR::EULER::
cal_repulsion_force_R_boost_with_RF
(POTENTIAL_SET& POTs, MATRIX& given_vec, long index_particle, RDIST& R_boost,
 double& RF_repulsion_xx, double& RF_repulsion_yy, double& RF_repulsion_zz,
 double& RF_repulsion_xy, double& RF_repulsion_xz, double& RF_repulsion_yz,
 double& energy_repulsive_potential)
{
  double time_st = omp_get_wtime();
  given_vec.set_value(0.);
  long cell_index_particle = R_boost.cell_index[index_particle];
  for(long k=0; k<R_boost.N_neighbor_cells; k++)
    {
      long cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(long p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          long index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.Rsca[index_particle](index_target);
          if (index_target != index_particle)
            {
              // H-D08/F05: force, virial and potential describe the same
              // state-dependent pair coupling. Scaling force alone breaks
              // -grad(U)=F; the unequal-population finite-difference fixture
              // detects a mismatch hidden by unit/pre-averaged coupling.
              const double coupling = POTs.f_repulsion_coupling(POTs.force_variables, index_particle, index_target);
              double repulsion = coupling*POTs.f_repulsion(distance, POTs.force_variables);
              // double repulsion = POTs.f_repulsion(distance, POTs.force_variables);
	      
              double repulsion_div_distance = repulsion/distance;
              MATRIX& rel_vector = R_boost.Rvec[index_particle][index_target];
              // for(long d=0; d<R_boost.N_dimension; d++)
              //   {
              //     given_vec(d) += repulsion*R_boost.Rvec[index_particle][index_target](d);
              //   }
              // the following make overhead.
              cblas_daxpy(given_vec.size,
                          repulsion/distance,
                          // R_boost.Rvec[index_particle][index_target].data,
                          rel_vector.data,
                          1,
                          given_vec.data,
                          1);

              RF_repulsion_xx += repulsion_div_distance*rel_vector(0)*rel_vector(0);
              RF_repulsion_yy += repulsion_div_distance*rel_vector(1)*rel_vector(1);
              RF_repulsion_zz += repulsion_div_distance*rel_vector(2)*rel_vector(2);

              RF_repulsion_xy += repulsion_div_distance*rel_vector(0)*rel_vector(1);
              RF_repulsion_xz += repulsion_div_distance*rel_vector(0)*rel_vector(2);
              RF_repulsion_yz += repulsion_div_distance*rel_vector(1)*rel_vector(2);

              energy_repulsive_potential += coupling*POTs.e_repulsion(distance, POTs.force_variables);
	      
            }
        }
    }
  return omp_get_wtime() - time_st;
}



double
INTEGRATOR::EULER::
cal_random_force_boost
(POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost)
{
  double time_st = omp_get_wtime();
  RANDOM::single_random_vector_generator_variance_boost(given_vec, 1.0, r_boost);
  matrix_mul(given_vec, sqrt(2.0));
  POTs.scale_random(given_vec, POTs.force_variables);
  return omp_get_wtime() - time_st;
}

double
INTEGRATOR::EULER::
cal_random_force_boost_simplified
(POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost)
{
  // this is duplicate one with INTEGRATOR::EULER::cal_random_force_boost
  // however, it simplified and reduce the overhead
  double time_st = omp_get_wtime();
  double variance = 1.0;
  double pre_factor = sqrt(3.0)*sqrt(2.0);
  for(long k=0; k<given_vec.size; k++)
    {
      given_vec(k) = pre_factor*variance*(2.0*(gsl_rng_uniform(r_boost) - 0.5));
    }
  return omp_get_wtime() - time_st;
}

double
ANALYSIS::
cal_total_energy_R_boost(POTENTIAL_SET& POTs, RDIST& R_boost)
{
  return cal_potential_energy_R_boost(POTs, R_boost);// + cal_kinetic_energy(TRAJ, POTs, index_t);
}

double
ANALYSIS::ANAL_ASSOCIATION::
cal_potential_energy_R_boost
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, RDIST& R_boost)
{
  double energy = 0.;

  for(long i=0; i<CONNECT.Np; i++)
    {
      for(long j=0; j<CONNECT.TOKEN[i]; j++)
        {
          energy += CONNECT.weight[i](j)*POTs.e_connector(R_boost.Rsca[i]((long)CONNECT.HASH[i](j)), POTs.force_variables);
        }
    }
  return energy/2. + ANALYSIS::cal_potential_energy_R_boost(POTs, R_boost);
}

double
ANALYSIS::DUMBBELL::
cal_potential_energy_R_boost
(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, RDIST& R_boost)
{
  double energy = 0.;
  for (long i=0; i<CONNECT.Np; i++)
    {
      for(long j=0; j<CONNECT.TOKEN[i]; j++)
	{
	  energy += POTs.e_connector(R_boost.Rsca[i]((long)CONNECT.HASH[i](j)), POTs.force_variables);
	}
    }
  return energy/2.; // each undirected dumbbell is visited from both endpoints
}

double
ANALYSIS::
cal_potential_energy_R_boost
(POTENTIAL_SET& POTs, RDIST& R_boost)
{
  double energy = 0.;
  for(long index_particle=0; index_particle<R_boost.Np; index_particle++)
    {
      long cell_index_particle = R_boost.cell_index[index_particle];
      for(long k=0; k<R_boost.N_neighbor_cells; k++)
        {
          long cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
          for(long p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
            {
              long index_target = R_boost(cell_index_neighbor, p);
              double distance = R_boost.Rsca[index_particle](index_target);
              if(index_target!=index_particle)
                energy += .5*POTs.f_repulsion_coupling(POTs.force_variables,index_particle,index_target)
                            *POTs.e_repulsion(distance,POTs.force_variables);
            }
        }
    }
  return energy;
}
double
ANALYSIS::
CAL_ENERGY_BROWNIAN
(POTENTIAL_SET& POTs, MATRIX& mat_energy, double time)
{
  double time_st = omp_get_wtime();
  mat_energy(0) = time;
  mat_energy(2) = 0.; // this is temporal setting
  // mat_energy(2) = cal_potential_energy_R_boost(POTs, R_boost);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return omp_get_wtime() - time_st;
}


double
ANALYSIS::
CAL_ENERGY_R_boost
(POTENTIAL_SET& POTs, MATRIX& mat_energy, double time, RDIST& R_boost)
{
  double time_st = omp_get_wtime();
  mat_energy(0) = time;
  mat_energy(2) = cal_potential_energy_R_boost(POTs, R_boost);
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return omp_get_wtime() - time_st;
}

double
ANALYSIS::DUMBBELL::
CAL_ENERGY_R_boost
(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& mat_energy, double time, RDIST& R_boost)
{
  double time_st = omp_get_wtime();
  mat_energy(0) = time;
  mat_energy(2) = ANALYSIS::DUMBBELL::cal_potential_energy_R_boost(POTs, CONNECT, R_boost);
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return omp_get_wtime() - time_st;
}

double
ANALYSIS::ANAL_ASSOCIATION::
CAL_ENERGY_R_boost
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, double time, RDIST& R_boost)
{
  /*
    Note that this code is old part of previous version.
    Now, the new analysis part is individually calculate during simulation because of efficiency
    The individual part will be found in the block of code for analysis in the main lines
   */
  // mat_energy(0) = (TRAJ.c_t - 1)*TRAJ.dt;
  double time_st = omp_get_wtime();
  mat_energy(0) = time;
  mat_energy(2) = ANALYSIS::ANAL_ASSOCIATION::cal_potential_energy_R_boost(POTs, CONNECT, R_boost);
  
  // the functionality for kinetic energy is disabled since the evolution equation on this code is using Weiner process that does not support differentiability for the position. On this regards, measuring the velocity cannot be obtained by this environment. For further detail, see the documents.
  mat_energy(1) = mat_energy(2) + mat_energy(3);
  return omp_get_wtime() - time_st;
}

