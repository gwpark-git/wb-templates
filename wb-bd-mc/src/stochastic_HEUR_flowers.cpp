#include "stochastic_HEUR_flowers.h"
#include <vector>

using namespace HEUR;

long
HEUR::
safety_random_micelle_selection(ASSOCIATION& CONNECT, gsl_rng* RNG_BOOST_SS_IT, long& Np)
{
  /*
    This safety function is design to avoid selection '0' chain end micelle as subject.
    The '0' chain end micelle is possible during simulation, but it is only acted by target micelle, i.e., suggestion connection junctions
    Note that even if this case, the adjacent list should be consistent with TOKEN 1, which means itself
   */
  long selected_micelle = -1;
  do
    {
      selected_micelle = RANDOM::return_LONG_INT_rand_boost(RNG_BOOST_SS_IT, Np);
      // if  (CONNECT.TOKEN[selected_micelle] == 1 && CONNECT.weight[selected_micelle](0) == 0)
      //   printf("IT HAPPENS!!\n");
    } while (CONNECT.TOKEN[selected_micelle] == 1 && CONNECT.weight[selected_micelle](0) == 0);
  return selected_micelle;
}

bool
HEUR::
get_index_bridge_distance_larger_than_cutoff_range(long& index_particle, long& index_hash, long& index_target, ASSOCIATION& CONNECT, RDIST& R_boost, double cutoff_range)
{
  
  for (long i=0; i<CONNECT.Np; i++)
    {
      for(long j=1; j<CONNECT.TOKEN[i]; j++)
	{
	  long tmp_index_target = CONNECT.HASH[i](j);
	  double distance = R_boost.Rsca[i](tmp_index_target);
	  if (distance > cutoff_range)//1.0)
	    {
	      index_particle = i;
	      index_hash = j;
	      index_target = tmp_index_target;
	      // printf("%d, %d, %d, %f\n", i, j, tmp_index_target, distance);
	      return 1;
	    }
	}
    }
  return 0;
}

double
HEUR::
record_RDIST_PARTICLE(ofstream &file_DIST,
                      RDIST& R_boost)
{
  double time_st = omp_get_wtime();

  for(long index_particle=0; index_particle<R_boost.Np/2; index_particle++)
    {
      long cell_index_particle = R_boost.cell_index[index_particle];
      for(long k=0; k<R_boost.N_neighbor_cells; k++)
        {
          long cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
          for(long p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
            {
              long index_target = R_boost(cell_index_neighbor, p);
              double distance = R_boost.Rsca[index_particle](index_target);
              // if(distance <= R_boost.cell_length && index_particle != index_target)
              //   file_DIST << std::scientific << distance << '\t';
              bool check_asymmetric_cells = 0;
              for(long n=0; n<R_boost.Rvec[index_particle][index_target].size; n++)
                check_asymmetric_cells = check_asymmetric_cells | (R_boost.Rvec[index_particle][index_target](n) <= R_boost.cell_length[n]);
              if (index_particle != index_target && check_asymmetric_cells)
                file_DIST << std::scientific << distance << '\t';
            }
	  file_DIST << std::endl;
        }
    }
  // file_DIST << '\n';
  return omp_get_wtime() - time_st;
}

double
HEUR::
record_RDIST_BRIDGE(ofstream &file_DIST,
                    RDIST& R_boost, ASSOCIATION& CONNECT)
{
  double time_st = omp_get_wtime();
  for(long i=0; i<CONNECT.Np; i++)
    {
      for(long j=1; j<CONNECT.TOKEN[i]; j++)
        {
          double distance = R_boost.Rsca[i](CONNECT.HASH[i](j));
          file_DIST << std::scientific << distance << '\t';
        }
    }
  // file_DIST << '\n';
  return omp_get_wtime() - time_st;
}

long
HEUR::
stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition)
{

  using namespace std;
  double time_st_simulation = omp_get_wtime();
  // printf("N_min : %d\n", CONNECT.N_min);
  printf("### STOCHASTIC SIMULATION FOR HEUR FLOWER ###\n");
  printf("SETTING simulation environment ...");
  TEMPORAL_VARIABLE_HEUR VAR(given_condition, TRAJ.rows);
  LOCK LOCKER(TRAJ.Np);
  RDIST R_boost(given_condition);
  
  long &cnt_cancel = VAR.cnt_arr[INDEX_MC::CANCEL], &cnt_add = VAR.cnt_arr[INDEX_MC::ADD], &cnt_del = VAR.cnt_arr[INDEX_MC::OPP_DEL], &cnt_mov = VAR.cnt_arr[INDEX_MC::MOV], &cnt_lock = VAR.cnt_arr[INDEX_MC::LOCK], &cnt_SS = VAR.cnt_arr[5];
  cnt_add = CONNECT.N_ASSOCIATION;
  
  if(given_condition("MC_LOG") == "TRUE")
    {
      DATA.MC_LOG << "00_cnt"<< '\t'  << "01_index_itself"<< '\t'  << "02_roll_dCDF"<< '\t'  << "03_hash_index_target"<< '\t'  << "04_index_target"<< '\t'  << "05_roll_dCDF_U"<< '\t'  << "06_index_k_new_target"<< '\t'  << "07_index_new_target"<< '\t'  << "08_TOKEN(i_NT)"<< '\t'  << "09_N_CHAIN_ENDS" << '\t'<< "10_N_CHAIN_ITSELF" << '\t' << "11_N_TOTAL_ASSOCIATION*2" << '\t' << "12_cnt_add"<< '\t'  << "13_cnt_mov"<< '\t'  << "14_cnt_del"<< '\t'  << "15_cnt_cancel" << '\t' << "16_cnt_lock" << endl;
    }
  printf("SET SIMULATION PARAMETERS ...");

  MATRIX energy(1, VAR.N_components_energy, 0.);
  // // 0: time step to write 1-3: energy, 4: NAS, 5: real time
  // // 6: (xx)[RF], 7: (yy)[RF], 8: (zz)[RF], 9: (xy)[RF], 10: (xz)[RF], 11:(yz)[RF]

  VAR.time_DIST +=
    REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, 0, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);
  if(R_boost.CELL_LIST_BOOST)
    for(long i=0; i<TRAJ.Np; ++i)
      GEOMETRY::compute_bonded_RDIST_particle(R_boost,CONNECT,i,TRAJ,0);
  
  // VAR.time_AN +=
  //   ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);

  RNG_BOOST RNG(given_condition);
  
  INDEX_MC *IDX_ARR = new INDEX_MC [VAR.N_THREADS_SS]; // it will call default constructor

  printf("DONE\nSTART SIMULATION\n\n");
  if(VAR.STEP_SHEAR)
    {
      long time_init = 0;
      
      GEOMETRY::
        set_box_shift_factor(VAR.shear_PBC_shift, R_boost,
                             VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis],
                             TRAJ.box_dimension[VAR.shear_axis]);
      GEOMETRY::
        apply_step_shear(TRAJ, time_init,
                         VAR.shear_axis, VAR.shear_grad_axis,
                         VAR.gamma_0, TRAJ.box_dimension[VAR.shear_grad_axis]);
    }

      // VAR.shear_PBC_shift = fmod(VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis], TRAJ.box_dimension[VAR.shear_axis]);
      // // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
      // R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
      // long central_standard = (long)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
      // R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;

      // GEOMETRY::
      //   set_box_shift_factor(double &shear_PBC_shift,
      //                        RDIST& R_boost,
      //                        double maximum_displacement,
      //                        double box_dimension_shear_axis, double box_dimension_shear_grad_axis)
      //   {
      //     /* // Original Form
      //     VAR.shear_PBC_shift = fmod(VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis], TRAJ.box_dimension[VAR.shear_axis]);
      //     // VAR.shear_PBC_shift = VAR.gamma_0*TRAJ.box_dimension[VAR.shear_grad_axis];
      //     R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
      //     long central_standard = (long)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
      //     R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
      //     */
          
      //     shear_PBC_shift = fmod(maximum_displacement, box_dimension_shear_axis);
      //     R_boost.map_to_central_box_image = fmod(shear_PBC_shift, box_dimension_shear_axis);
      //     long central_standard = (long)(2*R_boost.map_to_central_box_image/box_dimension_shear_axis);
      //     R_boost.map_to_central_box_image -= box_dimension_shear_axis*(double)central_standard;
      //   }

  
  if(VAR.SIMPLE_SHEAR)
    {
      printf("APPLY SIMPLE_SHEAR\n");
    }
  for(long t = 0; t<VAR.Nt-1; t++)
    {
      long index_t_now = t % VAR.N_basic;
      long index_t_next = (t+1) % VAR.N_basic;
      TRAJ(index_t_next) = TRAJ(index_t_now) + TRAJ.dt; // it will inheritance time step from the previous input file
      ++TRAJ.c_t;
      long cnt = 1;
      VAR.virial_initial();

      // printf("test1\n");
      if(VAR.SIMPLE_SHEAR)
        {
          double time_div_tau_R = TRAJ(index_t_now);
          /*
          // check the basic shift factor. Note that relative distance in the upper and lower box in shear gradient direction becomes Wi_R*L_D*t
          VAR.shear_PBC_shift = fmod(VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_R, TRAJ.box_dimension[VAR.shear_axis]);
          // the modulo scheme will works for efficiency of getting the directly connected boundary of PBC box
          R_boost.map_to_central_box_image = fmod(VAR.shear_PBC_shift, TRAJ.box_dimension[VAR.shear_axis]);
          
          // central_standard to identify the middle-point shift factor since something beyond middle point represent the longer upper (or lower) boxes.
          // it is of importance to identify minimum image convention
          long central_standard = (long)(2*R_boost.map_to_central_box_image/TRAJ.box_dimension[VAR.shear_axis]);
          // this is the last parts. The corrected re-shift the relative box distance based on middle point identification
          R_boost.map_to_central_box_image -= TRAJ.box_dimension[VAR.shear_axis]*(double)central_standard;
          // note that map_to_central_box_image will be used when we try to measure minimum distance between two micelles
          */
          GEOMETRY::
            set_box_shift_factor(VAR.shear_PBC_shift,
                                 R_boost,
                                 VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*time_div_tau_R,
                                 TRAJ.box_dimension[VAR.shear_axis]);
        }
      // printf("test1\n");      
      VAR.time_DIST +=
        REPULSIVE_BROWNIAN::OMP_compute_RDIST(TRAJ, index_t_now, R_boost, VAR.tmp_index_vec, VAR.N_THREADS_BD);
      // Bond geometry is independent of repulsive neighbor coverage.
      if(R_boost.CELL_LIST_BOOST) {
        double bonded_start=omp_get_wtime();
#pragma omp parallel for num_threads(VAR.N_THREADS_BD) if(VAR.N_THREADS_BD > 1)
        for(long i=0; i<TRAJ.Np; ++i)
          GEOMETRY::compute_bonded_RDIST_particle(R_boost,CONNECT,i,TRAJ,index_t_now);
        VAR.time_DIST+=omp_get_wtime()-bonded_start;
      }
      // printf("test2\n");      
      VAR.time_SS +=
        OMP_SS_topological_time_evolution(t, CONNECT, CHAIN, POTs, R_boost, RNG, IDX_ARR, DATA, given_condition, LOCKER, VAR);
      // printf("test3\n");      
      VAR.time_LV +=
        OMP_time_evolution_Euler(TRAJ, index_t_now, index_t_next,
                                 CONNECT,
                                 POTs, VAR.force_repulsion, VAR.force_random, VAR.force_spring,
                                 R_boost, VAR.vec_boost_Nd_parallel,
                                 RNG,
                                 VAR.N_THREADS_BD,
                                 given_condition, VAR);

      // printf("\n\n=====CHECK number of chain ends\n\n");
      // for(long i=0; i<TRAJ.Np; i++)
      // 	{
      // 	  if((long)CONNECT.N_ASSOCIATION_PARTICLE(i) != (long)CONNECT.get_N_association_particle(i))
      // 	    {
      // 	      printf("ERR N_ASSOCIATION_PARTICLE: p(%ld)=%ld (expected %ld)\n", i, (long)CONNECT.N_ASSOCIATION_PARTICLE(i), CONNECT.get_N_association_particle(i));
      // 	      for(long j=0; j<CONNECT.TOKEN[i]; j++)
      // 		printf("%ld\t", (long)CONNECT.weight[i](j));
      // 	      printf("\n");
      // 	    }
		       
      // 	}
      
      // printf("test4\n");
      if(VAR.SIMPLE_SHEAR || VAR.STEP_SHEAR)
        {
          double boundary_shift=VAR.shear_PBC_shift;
          if(VAR.SIMPLE_SHEAR)
            GEOMETRY::set_box_shift_factor(boundary_shift,
              VAR.Wi_tau_R*TRAJ.box_dimension[VAR.shear_grad_axis]*TRAJ(index_t_next),
              TRAJ.box_dimension[VAR.shear_axis]);
          // when the particles beyond boundary of PBC box in shear gradient direction, the minimum image convention should be changed
          // this functionality will remap in this case
          VAR.time_LV += 
            GEOMETRY::apply_shear_boundary_condition(TRAJ, index_t_next, VAR.shear_axis, VAR.shear_grad_axis, boundary_shift);
        }
      VAR.time_LV +=
        GEOMETRY::minimum_image_convention_loop(TRAJ, index_t_next); // applying minimum image convention for PBC
      // note that minimum_image_convention must be with loop when shear flow is implemented
      // this is due to the fact that the previous minimum_image_convention cannot cope the shift higher than two times of box dimension

      VAR.N_associations = cnt_add - cnt_del;

      if(t%VAR.N_skip_ener==0 || t%VAR.N_skip_file==0)
        {
          // VAR.time_AN +=
          //   ANALYSIS::ANAL_ASSOCIATION::CAL_ENERGY_R_boost(POTs, CONNECT, energy, (TRAJ.c_t - 1.)*TRAJ.dt, R_boost);
          energy(0) = TRAJ(index_t_now);
	  
          VAR.time_AN +=
            VAR.record_virial_into_energy_array(energy);
	  
          VAR.time_AN +=
            HEUR::sum_virial_components(energy);

          
          energy(4) = (double)VAR.N_associations; // number of associations
          energy(5) = omp_get_wtime() - time_st_simulation; // record computation time
          energy(30) = VAR.N_diff_associations; // it record number of different associations
          
          VAR.time_file +=
            energy.fprint_row(DATA.ener, 0);

          if(t%VAR.N_skip_file==0)
            {
              VAR.time_file +=
                record_simulation_data(DATA, TRAJ, CONNECT, CHAIN, index_t_now); // neeed review

              VAR.simulation_time = TRAJ(index_t_now)/(atof(given_condition("repulsion_coefficient").c_str())*atof(given_condition("Rt").c_str()));
          
              VAR.time_RECORDED +=
                report_simulation_info(TRAJ, energy, VAR);
              // H-D08/F12: disabled recording must not evaluate modulo an
              // uninitialized/zero cadence or write through a closed stream.
              // This guard is separate from the physical update schedule.
              if(DATA.r_dist_bridge.is_open() && VAR.N_skip_rdist > 0 && t%VAR.N_skip_rdist==0)
                {
                  VAR.time_RECORDED +=
                    record_RDIST_BRIDGE(DATA.r_dist_bridge, R_boost, CONNECT);
                  VAR.time_RECORDED +=
                    record_RDIST_PARTICLE(DATA.r_dist_particle, R_boost);
                }
            }

        }
    }

  double time_simulation = omp_get_wtime() - time_st_simulation;
  printf("Total simulation time = %6.3e, recorded time = %6.3e (record/total = %3.1f)\n", time_simulation, VAR.time_RECORDED, VAR.time_RECORDED/time_simulation);
  delete[] IDX_ARR;
  return 0;
}

double
HEUR::
record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, const long index_t_now)
{
  double time_st = omp_get_wtime();
  TRAJ.fprint_row(DATA.traj, index_t_now);
  
  for(long ip=0; ip<TRAJ.Np; ip++)
    {
      CONNECT.HASH[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.hash, 1, CONNECT.TOKEN[ip]);
      CONNECT.weight[ip].fprint_LONG_skip_transpose_LIMROWS(DATA.weight, 1, CONNECT.TOKEN[ip]);
    }

  if(CHAIN.INITIALIZATION)
    CHAIN.write(DATA.chain);
  return omp_get_wtime() - time_st;
}

double
HEUR::
report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE_HEUR& VAR)
{
  const double report_started=omp_get_wtime();
  printf(" PBC VOLUME = %3.2e\n", VAR.volume_PBC_box);
  double total_time = VAR.time_SS + VAR.time_LV + VAR.time_AN + VAR.time_file + VAR.time_DIST;
  double sum_time_LV = VAR.time_LV_init + VAR.time_LV_force + VAR.time_LV_update;
  printf("##### STEPS = %ld\tTIME = %8.6e tau_0\tENERGY = %6.3e (computing time: %4.3e)\n", TRAJ.c_t, VAR.simulation_time, energy(1), energy(5));
  printf("time consuming: LV = %3.2e (%3.1f%%), SS = %3.2e (%3.1f%%), AN = %3.2e (%3.1f%%), FILE = %3.2e (%3.1f%%), DIST = %3.2e (%3.1f%%)\n", VAR.time_LV, VAR.time_LV*100/total_time, VAR.time_SS, VAR.time_SS*100/total_time, VAR.time_AN, VAR.time_AN*100/total_time, VAR.time_file, VAR.time_file*100/total_time, VAR.time_DIST, VAR.time_DIST*100/total_time);
  printf("LV: init = %3.1f%%, force = %3.1f%% (repulsion = %3.1f%%, random = %3.1f%%, connector = %3.1f%%), update = %3.1f%%\n", 100.*VAR.time_LV_init/sum_time_LV, 100.*VAR.time_LV_force/sum_time_LV, 100.*VAR.time_LV_force_repulsion/VAR.time_LV_force, 100.*VAR.time_LV_force_random/VAR.time_LV_force, 100.*VAR.time_LV_force_connector/VAR.time_LV_force, 100.*VAR.time_LV_update/sum_time_LV);
  printf("SS: index process = %3.1f%%, LOCKING = %3.1f%%, check_dissociation = %3.1f%%, transition = %3.1f%%, update info. = %3.1f%%, pre-ASSOCIATION = %3.1f%%, pre-SUGGESTION = %3.1f%%)  \n", 100.*VAR.time_SS_index/VAR.time_SS, 100.*VAR.time_SS_LOCK/VAR.time_SS, 100.*VAR.time_SS_check/VAR.time_SS, 100.*VAR.time_SS_transition/VAR.time_SS, 100.*VAR.time_SS_update_info/VAR.time_SS, 100.*VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS, 100.*VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS);
  // 	 100.*VAR.time_SS_update/VAR.time_SS, VAR.time_SS_update_ASSOCIATION_MAP, VAR.time_SS_update_ASSOCIATION_MAP/VAR.time_SS_update, VAR.time_SS_update_CHAIN_SUGGESTION, VAR.time_SS_update_CHAIN_SUGGESTION/VAR.time_SS_update);
  printf("CHECK LAST STATISTICS: N_tot_asso = %ld, NAS = %ld, fraction=%4.3f, dNAS = %ld ####\n\n", VAR.N_tot_associable_chain, VAR.N_associations, VAR.N_associations/(double)VAR.N_tot_associable_chain, VAR.N_diff_associations);
  return omp_get_wtime()-report_started;
}


double
HEUR::
OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const long index_t_now, const long index_t_next,
                         ASSOCIATION& CONNECT,
                         POTENTIAL_SET& POTs, MATRIX* force_repulsion, MATRIX* force_random, MATRIX* force_spring,
                         RDIST& R_boost, MATRIX* vec_boost_Nd_parallel,
                         RNG_BOOST& RNG,
                         const long N_THREADS_BD,
                         COND& given_condition, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = omp_get_wtime();
  double time_LV_init = 0., time_LV_force = 0., time_LV_update = 0.;
  double time_LV_force_repulsion = 0., time_LV_force_random = 0., time_LV_force_connector = 0.;

  double RF_random_xx = 0., RF_random_yy = 0., RF_random_zz = 0.;
  double RF_random_xy = 0., RF_random_xz = 0., RF_random_yz = 0.;
  double RF_repulsion_xx = 0., RF_repulsion_yy = 0., RF_repulsion_zz = 0.;
  double RF_repulsion_xy = 0., RF_repulsion_xz = 0., RF_repulsion_yz = 0.;
  double RF_connector_xx = 0., RF_connector_yy = 0., RF_connector_zz = 0.;
  double RF_connector_xy = 0., RF_connector_xz = 0., RF_connector_yz = 0.;

  double energy_repulsive_potential = 0., energy_elastic_potential = 0.;

  long N_diff_associations = 0;
  // H-D08/F02: friction belongs to a particle's current bridge degree, not
  // its coordinate index. Freeze one validated value for every coordinate in
  // this BD update; heterogeneous_particle_friction uses zero-noise drift.
  std::vector<double> friction(TRAJ.Np);
  for(long i=0; i<TRAJ.Np; ++i) {
    friction[i]=POTs.zeta_particle(CONNECT.N_PARTICLE_ASSOCIATION(i),POTs.force_variables);
    if(!std::isfinite(friction[i]) || friction[i]<=0.)
      throw std::invalid_argument("particle friction must be positive and finite");
  }

#pragma omp parallel for default(none) if(N_THREADS_BD > 1)             \
  shared(TRAJ, index_t_now, index_t_next,                               \
         CONNECT, POTs, force_repulsion, force_random, force_spring,	\
         R_boost, vec_boost_Nd_parallel,  RNG, VAR, given_condition, friction)	\
  num_threads(N_THREADS_BD)                                             \
  reduction(+: time_LV_init, time_LV_force, time_LV_update, time_LV_force_repulsion, time_LV_force_random, time_LV_force_connector, \
            RF_random_xx, RF_random_yy, RF_random_zz,                   \
            RF_random_xy, RF_random_xz, RF_random_yz,                   \
            RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,          \
            RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,          \
            RF_connector_xx, RF_connector_yy, RF_connector_zz,          \
            RF_connector_xy, RF_connector_xz, RF_connector_yz,          \
            N_diff_associations, energy_repulsive_potential, energy_elastic_potential)
  
  for (long i=0; i<TRAJ.Np; i++)
    {
      long it = omp_get_thread_num(); // get thread number for shared array objects
      double time_st_init = omp_get_wtime();
      force_repulsion[i].set_value(0);
      force_random[i].set_value(0);
      force_spring[i].set_value(0);
      double time_st_force = omp_get_wtime();
      time_LV_init += time_st_force - time_st_init;
      time_LV_force_connector +=
        INTEGRATOR::EULER_ASSOCIATION::
        cal_connector_force_boost_with_RF(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca,
                                          RF_connector_xx, RF_connector_yy, RF_connector_zz,
                                          RF_connector_xy, RF_connector_xz, RF_connector_yz,
                                          energy_elastic_potential);
	
      // INTEGRATOR::EULER_ASSOCIATION::cal_connector_force_boost(POTs, CONNECT, force_spring[i], i, R_boost.Rvec, R_boost.Rsca);
	
      time_LV_force_repulsion +=
        INTEGRATOR::EULER::
        cal_repulsion_force_R_boost_with_RF(POTs, force_repulsion[i], i, R_boost,
                                            RF_repulsion_xx, RF_repulsion_yy, RF_repulsion_zz,
                                            RF_repulsion_xy, RF_repulsion_xz, RF_repulsion_yz,
                                            energy_repulsive_potential);
      // 	INTEGRATOR::EULER::cal_repulsion_force_R_boost(POTs, force_repulsion[i], i, R_boost);
      // time_LV_force_random +=
      // INTEGRATOR::EULER::cal_random_force_boost_simplified(POTs, force_random[i], RNG.BOOST_BD[it]);
      INTEGRATOR::EULER::cal_random_force_boost(POTs, force_random[i], RNG.BOOST_BD[it]);
      
      double time_st_update = omp_get_wtime();
      time_LV_force += time_st_update - time_st_force;
      
      for (long k=0; k<TRAJ.N_dimension; k++)
        {
	  const double zeta_corrector = friction[i];
	  // double zeta_corrector = POTs.zeta_particle(N_connectors_k, POTs.force_variables);
          // TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k)
          //   + (TRAJ.dt/zeta_corrector)*((1./POTs.force_variables[0])*force_spring[i](k) + force_repulsion[i](k))
          //   + sqrt(TRAJ.dt/zeta_corrector)*force_random[i](k);
          TRAJ(index_t_next, i, k) = TRAJ(index_t_now, i, k)
            + (TRAJ.dt/zeta_corrector)*((1./POTs.force_variables.repulsion_coefficient)*force_spring[i](k) + force_repulsion[i](k))
            + sqrt(TRAJ.dt/zeta_corrector)*force_random[i](k);

          
        }
      if(VAR.SIMPLE_SHEAR)
        TRAJ(index_t_next, i, VAR.shear_axis) += TRAJ.dt*VAR.Wi_tau_R*TRAJ(index_t_now, i, VAR.shear_grad_axis);

      RF_random_xx += TRAJ(index_t_now, i, 0)*force_random[i](0)/sqrt(TRAJ.dt);
      RF_random_yy += TRAJ(index_t_now, i, 1)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_zz += TRAJ(index_t_now, i, 2)*force_random[i](2)/sqrt(TRAJ.dt);

      RF_random_xy += TRAJ(index_t_now, i, 0)*force_random[i](1)/sqrt(TRAJ.dt);
      RF_random_xz += TRAJ(index_t_now, i, 0)*force_random[i](2)/sqrt(TRAJ.dt);
      RF_random_yz += TRAJ(index_t_now, i, 1)*force_random[i](2)/sqrt(TRAJ.dt);

      N_diff_associations += CONNECT.TOKEN[i] - 1; // 
        
      time_LV_update += omp_get_wtime() - time_st_update;
    }
  VAR.time_LV_init += time_LV_init;
  VAR.time_LV_force += time_LV_force;
  VAR.time_LV_update += time_LV_update;

  VAR.time_LV_force_repulsion += time_LV_force_repulsion;
  VAR.time_LV_force_random += time_LV_force_random;
  VAR.time_LV_force_connector += time_LV_force_connector;

  /*
    The following virial stress measurements are related with the time scale.
    In this simulation code, the repulsive time scale is used in order to non-dimensionalization of Langevin equation.
    Hence, the repulsive potential is computed without repulsion coefficient C and the other potential is affected by this repulsion coefficient C. 
    For instance, bridge interaction force is divided by C and the random force is divided by sqrt(C).
   */

  // allocationc omputed RF values into VAR
  VAR.RF_random_xx = RF_random_xx/VAR.volume_PBC_box; VAR.RF_random_yy = RF_random_yy/VAR.volume_PBC_box; VAR.RF_random_zz = RF_random_zz/VAR.volume_PBC_box;
  VAR.RF_random_xy = RF_random_xy/VAR.volume_PBC_box; VAR.RF_random_xz = RF_random_xz/VAR.volume_PBC_box; VAR.RF_random_yz = RF_random_yz/VAR.volume_PBC_box;

  // the following are changed in order tune the time scale in rightful way

  // VAR.RF_repulsion_xx = RF_repulsion_xx*POTs.force_variables[0]; VAR.RF_repulsion_yy = RF_repulsion_yy*POTs.force_variables[0]; VAR.RF_repulsion_zz = RF_repulsion_zz*POTs.force_variables[0];
  // VAR.RF_repulsion_xy = RF_repulsion_xy*POTs.force_variables[0]; VAR.RF_repulsion_xz = RF_repulsion_xz*POTs.force_variables[0]; VAR.RF_repulsion_yz = RF_repulsion_yz*POTs.force_variables[0];
  double duplication_divisor = 2.;
  VAR.RF_repulsion_xx = RF_repulsion_xx/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_yy = RF_repulsion_yy/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_zz = RF_repulsion_zz/(duplication_divisor*VAR.volume_PBC_box);
  VAR.RF_repulsion_xy = RF_repulsion_xy/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_xz = RF_repulsion_xz/(duplication_divisor*VAR.volume_PBC_box); VAR.RF_repulsion_yz = RF_repulsion_yz/(duplication_divisor*VAR.volume_PBC_box);

  VAR.RF_connector_xx = RF_connector_xx/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yy = RF_connector_yy/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_zz = RF_connector_zz/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box);
  VAR.RF_connector_xy = RF_connector_xy/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_xz = RF_connector_xz/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yz = RF_connector_yz/(POTs.force_variables.repulsion_coefficient*duplication_divisor*VAR.volume_PBC_box);

  
  VAR.energy_elastic_potential = energy_elastic_potential/(POTs.force_variables.repulsion_coefficient*duplication_divisor);

  
  // VAR.RF_connector_xx = RF_connector_xx/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yy = RF_connector_yy/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_zz = RF_connector_zz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box);
  // VAR.RF_connector_xy = RF_connector_xy/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_xz = RF_connector_xz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box); VAR.RF_connector_yz = RF_connector_yz/(POTs.force_variables[0]*duplication_divisor*VAR.volume_PBC_box);

  
  // VAR.energy_elastic_potential = energy_elastic_potential/(POTs.force_variables[0]*duplication_divisor);
  // the divisor C_rep (POTs.force_variables[0]) is used in order to tune time scale as tau_R

  VAR.energy_repulsive_potential = energy_repulsive_potential/duplication_divisor;

  VAR.N_diff_associations = (long)N_diff_associations/2; // divisor 2 for removing duplicate count
  
  return omp_get_wtime() - time_st;
}

double
HEUR::
write_MC_LOG_if_TRUE(bool flag_MC_LOG,
                     RECORD_DATA& DATA,
                     ASSOCIATION& CONNECT,
                     const INDEX_MC& IDX,
                     const long cnt, const long* cnt_arr,
                     const double rolling_dCDF, const double rolling_dCDF_U) 
{
  double time_st = omp_get_wtime();
  if(flag_MC_LOG)
    {
      long total_bonds = CONNECT.N_TOTAL_ASSOCIATION();
			  
      {
        DATA.MC_LOG << cnt << '\t' << IDX.beads[CONNECT.flag_itself] << '\t' << setprecision(7) << rolling_dCDF<< '\t'  << IDX.beads[CONNECT.flag_hash_other] << '\t'  << IDX.beads[CONNECT.flag_other] << '\t'  << setprecision(7) << rolling_dCDF_U<< '\t'  << IDX.beads[CONNECT.flag_hash_backtrace] << '\t'  << IDX.beads[CONNECT.flag_new] << '\t'  << CONNECT.TOKEN[IDX.beads[CONNECT.flag_itself]]<< '\t'<< CONNECT.N_CONNECTED_ENDS(IDX.beads[CONNECT.flag_itself]) << '\t' << CONNECT.weight[IDX.beads[CONNECT.flag_itself]](0) <<'\t' <<  total_bonds << '\t'  << cnt_arr[INDEX_MC::ADD]<< '\t'  << cnt_arr[INDEX_MC::MOV]<< '\t'  << cnt_arr[INDEX_MC::OPP_DEL]<< '\t'  << cnt_arr[INDEX_MC::CANCEL] << '\t' << cnt_arr[INDEX_MC::LOCK] << endl;
      }
    } // MC_LOG
  return omp_get_wtime() - time_st;
}

double
HEUR::
transition_single_chain_end(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, CHAIN_HANDLE& CHAIN, RDIST& R_boost, INDEX_MC& IDX, const long IDENTIFIER_ACTION)
{
  double time_st = omp_get_wtime();
  ACTION::ACT(POTs, CONNECT, IDX, R_boost.Rsca, IDENTIFIER_ACTION);
  if(CHAIN.INITIALIZATION)
    {
      // H-D07/H-D08: preserve the exact end after aggregate mutation. The
      // current serial attempt owns both representations and the subsequent
      // CDF/counter/log refresh; this helper alone is not an atomic transaction.
      CHAIN.TRACKING_ACTION(CONNECT, IDENTIFIER_ACTION, IDX); // it will track individual chain information
    }
  return omp_get_wtime() - time_st;
}

double
HEUR::
check_dissociation_probability(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, INDEX_MC& IDX, gsl_rng* RNG_BOOST_SS_IT, long& IDENTIFIER_ACTION)
{
  double time_st = omp_get_wtime();
  long index_itself = IDX.beads[CONNECT.flag_itself];
  long index_attached_bead = IDX.beads[CONNECT.flag_other];
  
  double distance_exist_bridge = R_boost.Rsca[index_itself](index_attached_bead);
  // double tpa = POTs.transition(distance_exist_bridge, POTs.f_connector(distance_exist_bridge, POTs.force_variables), POTs.force_variables);
  double tpa = POTs.transition(distance_exist_bridge, POTs.f_connector, POTs.force_variables);

  if (tpa == 1.0)
    {
      IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX);
      // if(IDX.beads[CONNECT.flag_itself] == 25 && IDX.beads[CONNECT.flag_new] == 25 && IDX.beads[CONNECT.flag_other] == 106)
      // 	{
      // 	  printf("tpa = %f, IDENTIFIER_ACTION = %d\n", tpa, IDENTIFIER_ACTION);
      // 	}      
    }
  else
    {
      double rolling_transition = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
      if (rolling_transition < tpa)
        {
          IDENTIFIER_ACTION = ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(CONNECT, IDX);
        }
      else
        IDENTIFIER_ACTION = IDX.CANCEL;
    }
  // if(IDX.beads[CONNECT.flag_itself] == 25 && IDX.beads[CONNECT.flag_new] == 25 && IDX.beads[CONNECT.flag_other] == 106)
  //   {
  //     printf("tpa = %f, IDENTIFIER_ACTION = %d\n", tpa, IDENTIFIER_ACTION);
  //   }
  return omp_get_wtime() - time_st;
}

double
HEUR::
LOCKING_PARALLEL(LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR, const INDEX_MC& IDX, long& IDENTIFIER_ACTION, long& IDENTIFIER_LOCKING)
{
  double time_st = omp_get_wtime();
  /*
    On the omp critical region, the block will work only one thread.
    If the other thread reaching this reason while there is one thread already working on this block, then the reached thread will wait until finishing the job of the other thread.
    This benefits to identify the working beads index on this case, since the 
  */
  // CHECKING

  for(long I_BEADS = 0; I_BEADS < 3; I_BEADS++)
    {
      if(LOCKER(IDX.beads[I_BEADS]))
        {
          IDENTIFIER_ACTION = IDX.CANCEL;
          IDENTIFIER_LOCKING = TRUE;
          break;
        }
    }

  // this is LOCKING procedure
  if(!IDENTIFIER_LOCKING)
    {
      VAR.cnt_SS ++;
      // cnt++;  // preventing LOCKING affect to the IDENTIFICATION of stochastic balance
      // for(long I_BEADS = 0; I_BEADS < 3 && N_THREADS_SS > 1; I_BEADS++) // N_THREADS_SS > 1 is not necessary
      for(long I_BEADS = 0; I_BEADS < 3; I_BEADS++)
        {
          LOCKER(IDX.beads[I_BEADS]) = TRUE;
        }
    }
  else
    {
      VAR.cnt_arr[INDEX_MC::LOCK] ++;
    }
  return omp_get_wtime() - time_st;
}

double
HEUR::
release_LOCKING(LOCK& LOCKER, INDEX_MC& IDX)
{
  double time_st = omp_get_wtime();
  for(long I_BEADS = 0; I_BEADS < 3; I_BEADS++)
    {
      LOCKER(IDX.beads[I_BEADS]) = FALSE;
    }
  return omp_get_wtime() - time_st;
}

double
HEUR::
micelle_selection(ASSOCIATION& CONNECT, gsl_rng* RNG_BOOST_SS_IT, INDEX_MC& IDX, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR, double& rolling_dCDF, double& rolling_dCDF_U, long& time_step)
{
  double time_st = omp_get_wtime();
  if (VAR.SIGN_DESTROY == TRUE)
    {
      long IDENTIFIER = get_index_bridge_distance_larger_than_cutoff_range(IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_hash_other], IDX.beads[CONNECT.flag_other], CONNECT, R_boost, VAR.CUTOFF_RANGE);
      if (IDENTIFIER == 0)
        {
          VAR.SIGN_DESTROY = FALSE;
        }
      else
        {
          time_step --; // it prevent to count the destroy process of bridge when it's distance is larger than 1.
        }
    }

  if (VAR.SIGN_DESTROY == FALSE)  // it must be checked in individual way since the first block can be skipped
    {
      // printf("tmp_check_block\n");
      // IDX.beads[CONNECT.flag_itself] = RANDOM::return_LONG_INT_rand_boost(RNG_BOOST_SS_IT, VAR.Np);
      IDX.beads[CONNECT.flag_itself] = safety_random_micelle_selection(CONNECT, RNG_BOOST_SS_IT, VAR.Np);
      // choice for selected chain end
      rolling_dCDF = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
      IDX.beads[CONNECT.flag_hash_other] = CONNECT.GET_INDEX_HASH_FROM_ROLL(IDX.beads[CONNECT.flag_itself], rolling_dCDF); 
      IDX.beads[CONNECT.flag_other] = CONNECT.HASH[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_hash_other]); 
      // choice for behaviour of selected chain end
    }
  // printf("1\n");
  rolling_dCDF_U = RANDOM::return_double_rand_SUP1_boost(RNG_BOOST_SS_IT);
  // printf("2\n");  
  // the PDF is already computed in the previous map
  IDX.beads[CONNECT.flag_hash_backtrace] = SEARCHING::backtrace_cell_list(CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]], CONNECT.TOKEN_ASSOCIATION[IDX.beads[CONNECT.flag_itself]], rolling_dCDF_U, IDX.beads[CONNECT.flag_itself], R_boost);
  // printf("3\n");
  // printf("%d, %d\n", IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_hash_backtrace]);
  IDX.beads[CONNECT.flag_new] = CONNECT.INDEX_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_hash_backtrace]);

  // if(IDX.beads[CONNECT.flag_itself] == 25 && IDX.beads[CONNECT.flag_new] == 25 && IDX.beads[CONNECT.flag_other] == 106)
  //   {
  //     double distance_other = R_boost.Rsca[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_other]);
  //     printf("current: %d, (%d, %d, %d, %d), %f, weight_itself=%d\n", time_step, IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_hash_other], IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_new], distance_other, (long)CONNECT.weight[IDX.beads[CONNECT.flag_itself]](0));
  //     // printf("dCDF: %f, %f\n", CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](0), CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](1));
  //     for(long i=0; i<CONNECT.TOKEN_ASSOCIATION[i]; i++)
  // 	{
  // 	  long index = CONNECT.Np - i - 1;
  // 	  printf("%d: %f\n", (long)CONNECT.INDEX_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](index), CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](index));
  // 	}
  //     printf("\n");

  //   }
  // printf("done\n");
  // double distance = R_boost.Rsca[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_new]);
  // if(distance > 1.0)
  //   {
  //     printf ("WRONG! %f, %f\n", distance, CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](IDX.beads[CONNECT.flag_hash_other]));
  //     // for(long i=CONNECT.Np - CONNECT.TOKEN_ASSOCIATION[CONNECT.flag_itself]; i<CONNECT.Np; i++)
  //     // 	{
  //     // 	  printf("(%f, %f)\t", R_boost.Rsca[IDX.beads[CONNECT.flag_itself]](CONNECT.dCDF_ASSOCIATION[IDX.beads[CONNECT.flag_itself]](i));
  //     // 	}
  //     // printf("\n");
      
  //   }
  return omp_get_wtime() - time_st;
}

double
HEUR::
OMP_SS_update_topology(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, CHAIN_HANDLE& CHAIN, RNG_BOOST& RNG, RECORD_DATA& DATA, INDEX_MC* IDX_ARR, LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR)
{
  double time_st = omp_get_wtime();
  double time_SS_index = 0., time_SS_LOCK = 0., time_SS_check = 0., time_SS_transition = 0., time_SS_update_info = 0.;
  // Mutable topology, CDFs, labelled-chain RNG/scratch and counters form one
  // transition. Use the serial reference order and stream 0 for every SS count.
  // Langevin forces and per-particle statistical maps remain parallel.
  // H-D08: endpoint locks alone did not protect proposal reads, tracking scratch
  // or logging. This conservative schedule intentionally changes old unsafe SS
  // ordering. Do not claim old parallel RNG replay or scalable MC from it.
  // Canceled/rejected attempts still consume their original attempt budget.
  (void)LOCKER;
  INDEX_MC& idx=IDX_ARR[0];
  for(long tp=0; tp<VAR.N_tot_associable_chain; ++tp) {
    double rolling_dCDF=0., rolling_dCDF_U=0.;
    time_SS_index+=micelle_selection(CONNECT,RNG.BOOST_SS[0],idx,R_boost,VAR,
                                    rolling_dCDF,rolling_dCDF_U,tp);
    long action=INDEX_MC::CANCEL;
    if(VAR.SIGN_DESTROY)
      action=ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST_NO_RESTRICTION(CONNECT,idx);
    else
      time_SS_check+=check_dissociation_probability(CONNECT,POTs,R_boost,idx,
                                                   RNG.BOOST_SS[0],action);
    time_SS_transition+=transition_single_chain_end(CONNECT,POTs,CHAIN,R_boost,idx,action);
    time_SS_update_info+=ACTION::UPDATE_INFORMATION(CONNECT,idx,VAR.cnt_arr,action);
    ++VAR.cnt_arr[action];
    ++VAR.cnt_SS;
    write_MC_LOG_if_TRUE(VAR.MC_LOG,DATA,CONNECT,idx,VAR.cnt_SS,VAR.cnt_arr,
                         rolling_dCDF,rolling_dCDF_U);
  }

  // adding the measured time to the variable structure
  VAR.time_SS_index += time_SS_index;
  VAR.time_SS_LOCK += time_SS_LOCK;
  VAR.time_SS_check += time_SS_check;
  VAR.time_SS_transition += time_SS_transition;
  VAR.time_SS_update_info += time_SS_update_info;

  return omp_get_wtime() - time_st;
}

double
HEUR::
OMP_SS_update_STATISTICS(ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN,
                         POTENTIAL_SET& POTs,
                         RDIST& R_boost,
                         TEMPORAL_VARIABLE_HEUR& VAR,
			 RECORD_DATA& DATA)
{
  double time_st = omp_get_wtime();
  long Np = VAR.Np;
  if (VAR.MC_renewal)
    {
      CONNECT.set_initial_condition();
      // Renewing the aggregate graph also renews labeled endpoints and reverse
      // maps. Keeping the previous chain state would leave two disagreeing
      // descriptions of one network (serialized_topology_and_renewal fixture).
      if(CHAIN.INITIALIZATION) CHAIN.allocate_existing_bridges(CONNECT);
      for(long i=0; i<4; i++)
        VAR.cnt_arr[i] = 0;
    }

  double time_update_ASSOCIATION_MAP = 0;
  double time_update_CHAIN_SUGGESTION = 0;
  
#pragma omp parallel for default(none) if(VAR.N_THREADS_BD > 1)         \
  shared(Np, POTs, CONNECT, R_boost)                                    \
  num_threads(VAR.N_THREADS_BD)                                         \
  reduction(+: time_update_ASSOCIATION_MAP, time_update_CHAIN_SUGGESTION)	
  
  
  for(long index_particle = 0; index_particle < Np; index_particle++)
    {
      // update suggestion map after time evolution of Langevin equation
      time_update_CHAIN_SUGGESTION +=
        CONNECT.update_CHAIN_SUGGESTION_MAP_particle(index_particle, POTs, R_boost);
      // update association map after time evolution of Langevin equation
      
      time_update_ASSOCIATION_MAP +=
        CONNECT.update_ASSOCIATION_MAP_particle(index_particle, POTs, R_boost);
    }

  if(VAR.MC_ASSOCIATION)
    {
      for(long index_particle = 0; index_particle < Np; index_particle++)
	{
	  long hash_init = Np - CONNECT.TOKEN_ASSOCIATION[index_particle] + 1;
	  for(long k=hash_init; k<Np; k++)
	    {
	      DATA.MC_ASSOCIATION << CONNECT.dCDF_ASSOCIATION[index_particle](k) << '\t';
	      DATA.RDIST_ASSOCIATION << R_boost.Rsca[index_particle](CONNECT.INDEX_ASSOCIATION[index_particle](k)) << '\t';
	      DATA.INDEX_ASSOCIATION << CONNECT.INDEX_ASSOCIATION[index_particle](k) << '\t';
	    }
	  DATA.MC_ASSOCIATION << std::endl;
	  DATA.RDIST_ASSOCIATION << std::endl;
	  DATA.INDEX_ASSOCIATION << std::endl;
	}
    }
  
  VAR.time_SS_update_ASSOCIATION_MAP += time_update_ASSOCIATION_MAP;
  VAR.time_SS_update_CHAIN_SUGGESTION += time_update_CHAIN_SUGGESTION;
  return omp_get_wtime() - time_st;
}


double
HEUR::
OMP_SS_topological_time_evolution(const long time_step_LV,
                                  ASSOCIATION& CONNECT,
                                  CHAIN_HANDLE& CHAIN,
                                  POTENTIAL_SET& POTs,
                                  RDIST& R_boost,
                                  RNG_BOOST& RNG,
                                  INDEX_MC* IDX_ARR,
                                  RECORD_DATA& DATA,
                                  COND& given_condition,
                                  LOCK& LOCKER,
                                  TEMPORAL_VARIABLE_HEUR& VAR)
// note that the TRJECTORY class dependency is deleted because of RDIST class.
{
  double time_st = omp_get_wtime();
  if(time_step_LV%VAR.N_steps_block == 0) // the equilibration functionality is disabled at this moment. (it will be seperated for future works)
    {
      // this is rearranged in order to use the previously updated information when MC_renewal is not turned on.
      // #pragma omp for
      // initialization_topological_update(CONNECT, POTs, given_condition, VAR);

      if(VAR.CUTOFF_BRIDGES == TRUE)
	{
	  VAR.SIGN_DESTROY = TRUE;
	}
      
      VAR.time_SS_update +=
        HEUR::OMP_SS_update_STATISTICS(CONNECT, CHAIN, POTs, R_boost, VAR, DATA);

      VAR.time_SS_CORE +=
        HEUR::OMP_SS_update_topology(CONNECT, POTs, R_boost, CHAIN, RNG, DATA, IDX_ARR, LOCKER, VAR);

      // long t1, t2, t3;
      // long IDENTIFIER = HEUR::get_index_bridge_distance_larger_than_unity(t1, t2, t3, CONNECT, R_boost);
      // if (IDENTIFIER > 0)
      //   printf("IDENTIFIER is NOT zero! %d, (%d, %d, %d), distance = %f\n", time_step_LV, t1, t2, t3, R_boost.Rsca[t1](t3));
      
      
    }
  return omp_get_wtime() - time_st;
}

HEUR::TEMPORAL_VARIABLE_HEUR::
TEMPORAL_VARIABLE_HEUR(COND& given_condition, long given_N_basic)
  : REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE(given_condition, given_N_basic)
{

  long N_dimension = atoi(given_condition("N_dimension").c_str());
  N_steps_block = atol(given_condition("N_steps_block").c_str());
  N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());

  time_SS = 0.;
  time_SS_CORE = 0.;
  time_SS_index = 0.;
  time_SS_LOCK = 0.;
  time_SS_check = 0.;
  time_SS_transition = 0.;
  time_SS_update_info = 0.;

  time_SS_update = 0.;
  time_SS_update_ASSOCIATION_MAP = 0.;
  time_SS_update_CHAIN_SUGGESTION = 0.;

  time_LV_force_repulsion = 0.;
  time_LV_force_random = 0.;
  time_LV_force_connector = 0.;
      
  for(long i=0; i<6; i++)
    cnt_arr[i] = 0;

  cnt_SS = 0;
      
  force_spring = new MATRIX [Np];
  for(long i=0; i<Np; i++)
    force_spring[i].initial(N_dimension, 1, 0.);

  N_tot_associable_chain = Np*atoi(given_condition("N_chains_per_particle").c_str());

  MC_renewal = FALSE;
  if(given_condition("MC_renewal")=="TRUE")
    MC_renewal = TRUE;

  MC_LOG = FALSE;
  if(given_condition("MC_LOG")=="TRUE")
    MC_LOG = TRUE;

  MC_ASSOCIATION = FALSE;
  if(given_condition("MC_ASSOCIATION_MAP") == "TRUE")
    MC_ASSOCIATION = TRUE;
  
  virial_initial();
  N_components_energy = 31; // last 30 index for N_diff_associations


  // for cutoff_scheme
  CUTOFF_BRIDGES = FALSE;
  SIGN_DESTROY = FALSE;
  if(given_condition("transition_probability") == "CUTOFF_DISSOCIATION")
    {
      CUTOFF_BRIDGES = TRUE;
      SIGN_DESTROY = TRUE;
      CUTOFF_RANGE = 1.0; // it is related with R0 cutoff dissociation rule
    }
  
  // if(given_condition("connector") == "Modified_Gaussian" || given_condition("connector") == "FENE")
  if(given_condition("connector") == "FENE") 
    {
      /*
        For Gaussian chain, the R_max value only affect to cutoff. It is design to prevent the equilibrium simulation when the bridge chain longer than half dimension of box, the length will be change during minimum image convention. It is of importance to aware that ratio_RM_R0 set with Ld/2 is not useful since it is already regarded as specific index set is replaced by some images. Therefore, it must have at least some tolerance factor such as Tf*Ld/2 with Tf<1.

        For finite extensible chain, this factor definitely related with maximally extendable length scale. Keep in mind that this factor is unit of R0. If scale_factor_chain (alpha) = 2.0, ratio_RM_R0 = 4 means maximally extendable chain is alpha*ratio_RM_R0 (8 in this example) times larger than end-to-end distance of chain. Therefore, cutoff_connection value must be set with ratio_RM_R0 value or lesser than it.

        Note that the cutoff_connection no matter what the equilibrium simulation or not. Even if the case for cell_list in equilibrium simulation, it is used the same parameter. However, keep rule the cutoff_connection value. It will destroy all the bridge beyond cutoff_connection which means we have to set lower value when CELL_LIST is applied and it should set lower than Ld/2 for nonequilibrium simulation (at this moment, nonequilibrium simulation with CELL_LIST is not applied).

        From now on, the cutoff_connection value must be positive all the non-equilibrium simulation since it is related with safety of statistics. Gaussian with this factor with higher Wi value might be suffer from it, which however, better than wrong distance contribution from bridges.

      */
      CUTOFF_BRIDGES = TRUE;
      SIGN_DESTROY = TRUE;
      CUTOFF_RANGE = atof(given_condition("cutoff_connection").c_str());
    }  

}


