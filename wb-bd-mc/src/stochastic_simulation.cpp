
#include <iostream>
#include "../lib/trajectory.h"		// TRAJECTORY class
#include "../lib/association.h"		// ASSOCIATION class
#include "../lib/handle_association.h"	// CHAIN class
#include "../lib/potential.h"		// POTENTIAL_SET class
#include "../lib/file_IO.h"		// RECORD_DATA class, COND class

#include <string>

#include "brownian.h"			// basic Brownian motion
#include "repulsive_brownian.h"		// EQUILIBRATION simulation
#include "stochastic_HEUR_flowers.h"	// stochastic simulation for HEUR flowers
#include "dumbbell_model.h"		// dumbbell model based on Brownian motion
using namespace std;

int
help()
{
  cout << "USAGE of Brownian Dynamics with Quadratic Repulsive Potential\n";
  cout << "argv[1] == condition files including variables:\n";
  cout << "\tLONG INTEGER format: N_dimension, Nt, Np\n";
  cout << "\tDOUBLE FLOAT format: dt, repulsion_coefficient, effective_distance\n";
  cout << "\tSTRING format: filename_trajectory, filename_energy, filename_energy_detail\n";
  return 0;
}


int
main
(int argc, char* argv[])
{
  try {
  /* main print */
  cout << "STOCHASTIC SIMULATION for ASSOCIATING TELECHELIC POLYMERS" << endl;
  cout << "Version: " << VERSION << endl;
  cout << "Code developer: Gun Woo Park (gunwoo.park@unina.it)" << endl;
    
  if(argc==1)
    {
      help();
      return 0;
    }
  else
    {
        
      // inp_cnt started with index 1 since 0 is itself
      // single input file case, argc==2, argv[0] == execution, argv[1] == inputfile
      for(long inp_cnt = 1; inp_cnt < argc; inp_cnt ++)
        {
          /*
            This for-loop is deisgn to perform sequence of jobs using the given input file.
            The sequence is exactly the same with the given index for arguments.
            Note that the all the controls should be properly set with the given input files, which means the inheritance from the previous computation in the given input file is exactly matched with the output file name from the previous input file.
          */
          
          COND given_condition(argv[inp_cnt]);
          given_condition.cond_print();
          if(atol(given_condition("N_dimension").c_str())!=3)
            throw std::invalid_argument("legacy simulator requires N_dimension=3 (six stress components)");
          if(atol(given_condition("Np").c_str())<=0 || atol(given_condition("Nt").c_str())<2 ||
             atol(given_condition("N_THREADS_BD").c_str())<=0 ||
             atol(given_condition("N_skip_ener").c_str())<=0 ||
             atol(given_condition("N_skip_file").c_str())<=0)
            throw std::invalid_argument("particle/time/thread counts and output cadence must be positive");
          if(given_condition("Method")=="NAPLE_ASSOCIATION" &&
             (atol(given_condition("N_THREADS_SS").c_str())<=0 ||
              atol(given_condition("N_steps_block").c_str())<=0 ||
              atol(given_condition("N_chains_per_particle").c_str())<=0))
            throw std::invalid_argument("network threads, block length and chain count must be positive");
          if(given_condition("Integrator")!="Euler")
            throw std::invalid_argument("legacy simulator supports Integrator=Euler");


          const bool brownian_clock=given_condition("Method")=="BROWNIAN" ||
                                     given_condition("Method")=="DUMBBELL";
          const double step=atof(given_condition(brownian_clock ? "dt/tauB" : "dt/tauR").c_str());
          if(!std::isfinite(step) || step<=0.)
            throw std::invalid_argument("simulation timestep must be positive and finite");
          for(long axis=0;axis<3;++axis) {
            const double length=HANDLE_COND::get_LBk(given_condition,axis);
            if(!std::isfinite(length) || length<=0.)
              throw std::invalid_argument("box lengths must be positive and finite");
          }
          const bool shear=given_condition("MECHANICAL_PERTURBATION")=="SIMPLE_SHEAR" ||
                           given_condition("MECHANICAL_PERTURBATION")=="STEP_SHEAR";
          if(shear && (atoi(given_condition("shear_axis").c_str())!=0 ||
                       atoi(given_condition("shear_grad_axis").c_str())!=1))
            throw std::invalid_argument("legacy shear supports shear_axis=0 and shear_grad_axis=1");

          if(shear) GEOMETRY::validate_shear_box(HANDLE_COND::get_LBk(given_condition,0),
                                                HANDLE_COND::get_LBk(given_condition,1),
                                                HANDLE_COND::get_LBk(given_condition,2));

          long N_basic = 2;
          if (given_condition("Integrator") == "Euler") 
            {
              // at the moment, the only method available is simple Euler integrator
              // for this reason, this condition will not be varied during condition test
              // however, the explicit description should be needed for condition data
              // that is because to prevent potential errors when we implement higher order methods.
              // Note that the only needed information for Euler is the previous time step. Therefore, we need previous and now which is in 2.
              N_basic = 2;

            }

          TRAJECTORY TRAJ(given_condition, N_basic);
          POTENTIAL_SET POTs;
          RECORD_DATA DATA(given_condition);
          if(given_condition("Method") == "NAPLE_ASSOCIATION")
            {
              // note that the time scale for Langevin dynamics on here is affected by repulsive coefficient, C_rep, which have differ non-dimensional form compared with pure Brownian motion
              // printf("settings\n");
              if(given_condition("Step") == "EQUILIBRATION")
                {
                  if(FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition)!=0)
                    throw std::invalid_argument("invalid potential or kinetics selection");
                  REPULSIVE_BROWNIAN::main_EQUILIBRATION(TRAJ, POTs, DATA, given_condition);
                }
              else
                {
                  ASSOCIATION CONNECT(given_condition);

                  if(FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition)!=0)
                    throw std::invalid_argument("invalid potential or kinetics selection");
                  CHAIN_HANDLE CHAIN(given_condition, CONNECT);
                  POTs.force_variables.N_ASSOCIATION_PARTICLE = &CONNECT.N_ASSOCIATION_PARTICLE;
                  /*
                   this is added on purpose
                   the initialization outside the class constructor is weird
                   the main reason to use non-recommendable initialization for force_variables is due to the hierachical library load.
                  */
                  HEUR::stochastic_simulation_HEUR_flowers(TRAJ, POTs, CONNECT, CHAIN, DATA, given_condition);
                  // }
                }
            }
          else if(given_condition("Method") == "REPULSIVE_BROWNIAN")
            {
              if(FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTs, given_condition)!=0)
                    throw std::invalid_argument("invalid potential or kinetics selection");
              REPULSIVE_BROWNIAN::main_EQUILIBRATION(TRAJ, POTs, DATA, given_condition);
            }
          else if(given_condition("Method") == "BROWNIAN")
            {
              
              // As mentioned in repulsive Brownian motion, here the time scale is used as the basic Brownian motion. Therefore, the conditional file has time step dt/tauB instead of dt/tauR.
              if(FORCE::BROWNIAN::MAP_potential_set(POTs, given_condition)!=0)
                    throw std::invalid_argument("invalid potential or kinetics selection"); // null allocator
              BROWNIAN::main_PURE_BROWNIAN(TRAJ, POTs, DATA, given_condition);
            }
          else if(given_condition("Method") == "DUMBBELL")
            {
              // CONNECTIVITY CONNECT(given_condition); // note that
              // HASH includes the self entry at slot 0 and the partner at slot 1.
              long N_max_connection = 2;
              printf("test:init\n");
              CONNECTIVITY CONNECT(TRAJ.Np, N_max_connection);
              printf("test:aft:connect\n");
              if (DUMBBELL::generate_dumbbell_connectivity(CONNECT) != 0)
                return 1;
              printf("test:aft:dumbbell_connect\n");
              if(FORCE::DUMBBELL::MAP_potential_set(POTs, given_condition)!=0)
                    throw std::invalid_argument("invalid potential or kinetics selection");
              printf("test:aft:potential_set\n");
              DUMBBELL::main_DUMBBELL(TRAJ, CONNECT, POTs, DATA, given_condition);
            }
          else
            {
              throw std::invalid_argument("unknown simulation Method: "+given_condition("Method"));
            }
        }
    }
  return 0;
  } catch(const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}



