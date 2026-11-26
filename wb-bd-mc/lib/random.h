
#ifndef RANDOM_H
#define RANDOM_H

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}

#include <math.h>
#include "numerics.h"
#include <iostream>
#include "matrix.h"
#include "file_IO.h"

// note that the RANDOM package should not has the dependency with lib_traj.h
class RNG_BOOST
{
 public:
  long N_THREADS_BD;
  long N_THREADS_SS;
  gsl_rng **BOOST_BD;
  gsl_rng **BOOST_SS;
  bool INITIALIZATION_BD;
  bool INITIALIZATION_SS;
  RNG_BOOST()
    {
      std::cout<< "ERR: basic constructor for RNG_BOOST class is not allowed\n";
    }
  RNG_BOOST
    (COND& given_condition)
    {
      INITIALIZATION_BD = TRUE;
      const gsl_rng_type *T_boost;
      N_THREADS_BD = atol(given_condition("N_THREADS_BD").c_str());
      BOOST_BD = new gsl_rng* [N_THREADS_BD];
      gsl_rng_env_setup();
      T_boost = gsl_rng_default;
      for(long i=0; i<N_THREADS_BD; i++)
        {
          BOOST_BD[i] = gsl_rng_alloc(T_boost);
          long seed_BD = atol(given_condition("basic_random_seed").c_str());
          gsl_rng_set(BOOST_BD[i], seed_BD + i); // it set the seed with index i
        }
      INITIALIZATION_SS = FALSE; // default
      if(given_condition("Method") == "NAPLE_ASSOCIATION")
	if(given_condition("Step") != "EQUILIBRATION")
	 /* given_condition("Step")!="EQUILIBRATION" && given_condition("Method")=="NAPLE_ASSOCIATION") */
        // identify the NAPLE_ASSOCIATION topological evolution with equilibration condtion or not
	  {
	    INITIALIZATION_SS = TRUE;
	    N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
	    BOOST_SS = new gsl_rng* [N_THREADS_SS];
	    for(long i=0; i<N_THREADS_SS; i++)
	      {
		BOOST_SS[i] = gsl_rng_alloc(T_boost);
		long seed_SS = atol(given_condition("basic_random_seed_SS").c_str());
		gsl_rng_set(BOOST_SS[i], seed_SS + i); 
	      }
        }
    }

  ~RNG_BOOST()
    {
      if(INITIALIZATION_BD)
        {
          for(long i=0; i<N_THREADS_BD; i++)
            {
              gsl_rng_free(BOOST_BD[i]);
            }
          delete[] BOOST_BD;
          if(INITIALIZATION_SS)
            {
              for(long i=0; i<N_THREADS_SS; i++)
                {
                  gsl_rng_free(BOOST_SS[i]);
                }
              delete[] BOOST_SS;
            }
        }
    }
};


namespace RANDOM
{
  long
    random_vector_generator
    (MATRIX& R_VEC_TRANS);
  long
    single_random_vector_generator
    (MATRIX& given_vec);
  long
    single_random_vector_generator_variance
    (MATRIX& given_vec, double s_2);
  long
    single_random_vector_generator_boost
    (MATRIX& given_vec, gsl_rng* r_boost);
  long
    single_random_vector_generator_variance_boost
    (MATRIX& given_vec, double s_2, gsl_rng* r_boost);
  long
    single_unit_random_vector_generator
    (MATRIX& given_vec);
  long
    unit_random_vector_generator
    (MATRIX& R_VEC_TRANS);
  long
    unit_random_vector_generator_2D
    (MATRIX& R_VEC_TRANS);

  long
    return_LONG_INT_rand
    (long SUP);
  long
    return_LONG_INT_rand_boost
    (gsl_rng* r, long SUP);
  double
    return_double_rand_SUP1();
  double
    return_double_rand_SUP1_boost
    (gsl_rng* r);
  long
    get_LONG_ARR_rand_boost
    (gsl_rng* r, long SUP, long* given_long_arr, long N_arr);
  long
    get_DOUBLE_ARR_rand_boost
    (gsl_rng* r, double* given_double_arr, long N_arr);
  
}


#endif 
