
// Trajectory class definition for first aim
// As test purpose, it include 2-dimensional case, then extend to the 3-dimensional geometry.
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include "matrix.h"
#include "random.h"
// #include <omp.h>

/* #include "read_file_condition.h" */
#include "file_IO.h"
class TRAJECTORY : public MATRIX
{

 public:
  enum {x = 0, y = 1, z = 2}; // this is for convenience

  long init_ref();
  enum {N_t_test = 10, N_p_test = 3};
  // as mentioned in previously, the dimensionality on this class is set to 2 for testing purpose. In order to make generalized code, it will changed and implemented into initiator of this class, later.

  string Method;
  /* long N_energy_frequency; */
  
  long Nt;
  long Np;
  // geometry. Note that all the proposed scheme should be defined in dimensionless scheme.
  long N_dimension = 0;
  double *box_dimension = nullptr;

  // time step
  long c_t;
  double dt;


  // member function
  long
    Rt_i
    (long t, long i_bead, long direction)
  {
    return index(t, i_bead*2*N_dimension + 1 + direction);
  }
  
  long
    Vt_i
    (long t, long i_bead, long direction)
  {
    return index(t, i_bead*2*N_dimension + 1 + N_dimension + direction);
  }
  
  
  // initiator
  long
    initialization
    (long N_time, long N_particle, double given_dt);
  long
    initialization_COND
    (COND& given_condition);

 TRAJECTORY();
 TRAJECTORY
   (const TRAJECTORY& TRAJ); // its for copy-constructor. It must be classified for openmp variable
 TRAJECTORY& operator=(const TRAJECTORY& TRAJ);
 TRAJECTORY
   (COND& given_condition, long N_basic);

 /* // inlined */
 /* // operator overloading */
 /* // simple operator for time */
 long
   read_exist_traj
 (const char* fn_given_traj, long N_steps);
 double&
   operator()
   (long time_t)
 {
   // it will return the reference of time
   return data[index(time_t, 0)];
 }

 double&
   operator()
   (long time_t, long bead_i, long dimension_k)
 {
   long index_position = 2*N_dimension*bead_i + 1 + dimension_k;
   return data[index(time_t, index_position)];
 }

 double&
   operator()
   (long i_RV, long time_t, long bead_i, long dimension_k)
 {
   long index_position = 2*N_dimension*bead_i + 1 + dimension_k + N_dimension*i_RV;
   return data[index(time_t, index_position)];
 }
 
 virtual
   ~TRAJECTORY()
    {
      delete[] box_dimension;
      // The MATRIX and REF_MATRIX objects are automatically killed-by virtual destructor option of the library. This fact is tested with packages.
    }

 // data loader
};



long
traj_count_line
(char fn[]);

namespace GENERATOR
{
  long
    random_position_generator
    (TRAJECTORY& TRAJ);
  long
    random_position_generator_REF
    (TRAJECTORY& TRAJ, MATRIX& R_VEC_TRANS);
}



class TRAJECTORY_HDF5
{
 public:
  double ***data = nullptr;
  long Np;
  long N_dimension;
  double *box_dimension;

  long c_t;
  double dt;
  
  bool INITIALIZATION = false;
  long N_basic;
  
  TRAJECTORY_HDF5()
    {
      std::cout << "ERR: TRAJECTORY_HDF5 class is not support basic constructor\n";
    }
  TRAJECTORY_HDF5(COND& given_condition, long given_N_basic)
    {
      Np = atoi(given_condition("Np").c_str());
      N_dimension = atoi(given_condition("N_dimension").c_str());
      N_basic = given_N_basic; // it report how many time step should be stored in memory
      c_t = 0;
      dt = atof(given_condition("dt/tauR").c_str());
      box_dimension = new double [N_dimension];
      for(long k=0; k<N_dimension; k++)
        {
          box_dimension[k] = atof(given_condition("box_dimension").c_str());
        }
      data = new double** [N_basic];
      for(long t=0; t<N_basic; t++)
        {
          data[t] = new double* [Np];
          for(long i=0; i<Np; i++)
            {
              data[t][i] = new double [N_dimension];
            }
        }
      INITIALIZATION = TRUE;
    }
  ~TRAJECTORY_HDF5()
    {
      if(INITIALIZATION)
        {
          for(long t=0; t<N_basic; t++)
            {
              for(long i=0; i<Np; i++)
                {
                  delete[] data[t][i];
                }
              delete[] data[t];
            }
          delete[] data;
          delete[] box_dimension;
          INITIALIZATION = FALSE;
        }
    }
  
};

#endif

