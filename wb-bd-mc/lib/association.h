
#ifndef ASSOCIATION_H
#define ASSOCIATION_H

#define TRUE 1
#define FALSE 0

#include <iostream>
// #include <omp.h>
#include "numerics.h"

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}
#include "matrix.h"
#include "connectivity.h"
#include "potential.h"
/* #include "trajectory.h" */
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "geometry.h"

class ASSOCIATION
: public CONNECTIVITY
{
 public:
  long Nc; // number of chains per micelle
  long Tec; // tolerance of end chains
  long N_min; // 2Nc - Tec
  long N_max; // 2Nc + Tec
  long N_ASSOCIATION;
  bool MULTIPLE_CONNECTIONS;
  bool NUMBER_RESTRICTION;
  bool PROBABILITY_AREA;

  MATRIX N_ASSOCIATION_PARTICLE;
  MATRIX *weight;
  // related with suggestion probability (selecting a degenerated pair)
  MATRIX *CASE_SUGGESTION;
  MATRIX *dPDF_SUGGESTION;
  MATRIX *dCDF_SUGGESTION;
  double *Z_SUGGESTION;

  // related with association probability
  MATRIX *dCDF_ASSOCIATION;
  /* size_t **INDEX_ASSOCIATION; */
  MATRIX *INDEX_ASSOCIATION;
  long *TOKEN_ASSOCIATION;

  // member variables
  
  /* long flag_other, flag_new, flag_itself, flag_hash_other; */
  static const long flag_other = 0;
  static const long flag_new = 1;
  static const long flag_itself = 2;
  static const long flag_hash_other = 3;
  static const long flag_hash_backtrace = 4;
  
  bool
    (*ADD_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    (*OPP_DEL_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    (*CANCEL_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    (*CHECK_N_ADD_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    (*CHECK_N_MOV_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    (*CHECK_N_OPP_DEL_ASSOCIATION)
    (ASSOCIATION& CONNECT, long index_set[]);
  
  long
    N_TOTAL_ASSOCIATION
    ();

  long
  N_PARTICLE_ASSOCIATION
  (long index_particle);
  
  long
    N_CONNECTED_ENDS
    (long given_index);
  long
    N_TOTAL_CONNECTED_ENDS
    ();
  bool
    CHECK_EXIST
    (long index_A, long index_B);
  bool
    CHECK_EXIST_1D
    (long index_A, long index_B);
  long
    add_association
    (const long index_particle, long index_target);
  long
    add_association_INFO
    (POTENTIAL_SET& POTs, const long index_particle, long index_target, double distance);
  long
    opp_del_association
    (const long index_particle, long index_target);
  long
    opp_del_association_hash
    (const long index_particle, long hash_index_target);
  long
    opp_del_association_IK
    (long index_I, long hash_index_K);
  long
    opp_del_association_grab_IK
    (long index_I, long hash_index_K);
  
  long
    GET_INDEX_HASH_FROM_ROLL
    (const long index_particle, double rolled_p);
  long
    GET_HASH_FROM_ROLL
    (const long index_particle, double rolled_p);
  long
    FIND_HASH_INDEX
    (long index_A, long index_B);
  long
  update_N_association
  ();
  
  long
  get_N_association_particle
  (long index_particle);

  // At the moment, there is no namespace dependence between KINETICS::NORMALIZED and KINETICS::METROPOLIS.
  // The dependency is set with the w_function and transtion_probability that is developed with lib_potential.
  // For better performance, it would be better to tune the following with specific dependencies (future works).
  double
    update_CASE_SUGGESTION_particle_hash_target
    (POTENTIAL_SET& POTs, const long index_particle, long hash_index_target, double distance);

  double
    update_CASE_SUGGESTION_particle_target
    (POTENTIAL_SET& POTs, const long index_particle, long index_target, double distance);

  double
    update_Z_SUGGESTION_particle
    (const long index_particle);
  double
    update_dPDF_SUGGESTION_particle
    (const long index_particle);
  double
    update_dCDF_SUGGESTION_particle
    (const long index_particle);
  double
    update_CHAIN_SUGGESTION_MAP_particle
    (const long index_particle, POTENTIAL_SET& POTs, RDIST& R_boost);

  double
    update_ASSOCIATION_MAP_particle
    (const long index_particle, POTENTIAL_SET& POTs, RDIST& R_boost);

  long
    return_multiplicity
    (const long index_particle, const long target_particle)
  // this additional function is useful for tracking individual chain
  {
    for(long j=0; j<TOKEN[index_particle]; j++)
      {
        if(HASH[index_particle](j) == target_particle)
          return weight[index_particle](j);
      }
    return 0;
  }
  
 ASSOCIATION()
   : CONNECTIVITY(){}

  ASSOCIATION
    (COND& given_condition);
  
  ASSOCIATION
    (long number_of_particles, long number_of_chains_per_particles, long tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS);

  long
    dynamic_alloc
    ();
  long
    set_initial_condition
    ();
  long
    initial
    ();
  long
    initial_inheritance
    ();
  virtual
    ~ASSOCIATION()
    {
      delete[] weight;
      
      // related with suggestion probability
      delete[] CASE_SUGGESTION;
      delete[] dPDF_SUGGESTION;
      delete[] dCDF_SUGGESTION;
      delete[] Z_SUGGESTION;
      // related with association probability
      delete[] dCDF_ASSOCIATION;
      delete[] INDEX_ASSOCIATION;
      
      delete[] TOKEN_ASSOCIATION;
    }

  long
    read_exist_weight
  (const char* fn_weight, long N_steps);
  
};

long
index_set_val
(size_t* index, long given_val, long size_of_arr);


namespace TRUTH_MAP
{
  namespace MULTIPLE
  {
    bool
      CHECK_N_ADD_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);
    bool
      CHECK_N_OPP_DEL_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);
    bool
      CHECK_N_MOV_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);

    namespace NO_RESTRICTION // this option is used when there are no restriction for number of attached chain ends
    {
      bool
	return_TRUE
	(ASSOCIATION& CONNECT, long index_set[]);

    }
  }
  
  namespace SINGLE
  {
    bool
      CHECK_N_ADD_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);
    bool
      CHECK_N_MOV_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);
    bool
      CHECK_N_OPP_DEL_BOOST
      (ASSOCIATION& CONNECT, long index_set[]);
  }
  bool
    CANCEL_ASSOCIATION_BOOST
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    ADD_ASSOCIATION_BOOST
    (ASSOCIATION& CONNECT, long index_set[]);
  bool
    OPP_DEL_ASSOCIATION_BOOST
    (ASSOCIATION& CONNECT, long index_set[]);
  
  bool
    OPP_DEL_ASSOCIATION_BASIC
    (ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  
  bool
    NEW_ASSOCIATION_BASIC
    (ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool
    NEW_ASSOCIATION_SINGLE
    (ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  
  bool
    MOV_ASSOCIATION_BASIC
    (ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool
    MOV_ASSOCIATION_SINGLE
    (ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
}

namespace MICELLE
{
  double
  update_SOFT_REPULSION_P2(POTENTIAL_SET& given_POT, ASSOCIATION& CONNECT);

  double
  update_SOFT_REPULSION_P2_particle(POTENTIAL_SET& given_POT, ASSOCIATION& CONNECT, long index_particle);

}

namespace KINETICS
{
  double
    CONNECTIVITY_update_dCDF_SUGGESTION_particle
    (ASSOCIATION* CONNECT, const long index_particle);
  double
    CONNECTIVITY_update_CASE_SUGGESTION_particle_hash_target
    (ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const long index_particle, long hash_index_target, double distance);
  double
    CONNECTIVITY_update_CASE_SUGGESTION_particle_target
    (ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const long index_particle, long index_target, double distance);

  double
    CONNECTIVITY_update_Z_SUGGESTION_particle
    (ASSOCIATION* CONNECT, const long index_particle);
  double
    CONNECTIVITY_update_dPDF_SUGGESTION_particle
    (ASSOCIATION* CONNECT, const long index_particle);
}


class CHAIN_INFORMATION
{
  /*
    The CHAIN class is using CHAIN_NODE without inheritance.
    It contains dynamically allocated CHAIN_NODE.
  */
 public:

  CHAIN_NODE *CHAIN = nullptr;
  long N_chains = 0;
  long N_particles = 0;
  
  bool INITIALIZATION = false;
  
  long&
    HEAD
    (long given_chain_index)
    {
      return CHAIN[given_chain_index].HEAD;
    }
  long&
    TAIL
    (long given_chain_index)
    {
      return CHAIN[given_chain_index].TAIL;
    }

  /*
    On the same reason, the ATTACHED will return pointer variable rather than reference variable. This break the existing interface.
  */
  long&
    ATTACHED
    (long given_chain_index, long flag_HEAD_TAIL)
    {
      return CHAIN[given_chain_index].index(flag_HEAD_TAIL);
    }
  
  
  long
    mov_attachment
    (long target_particle, long given_chain_index, long flag_HEAD_TAIL)
  {
    /*
      flag_HEAD_TAIL is 0 for HEAD while 1 for TAIL.
    */
    ATTACHED(given_chain_index, flag_HEAD_TAIL) = target_particle;

    return 0;
  }

  long
    add_attachment
    (const long index_particle, long given_chain_index, long flag_HEAD_TAIL)
  {
    // there is no distingushable between add and mov attachment for chain node
    mov_attachment(index_particle, given_chain_index, flag_HEAD_TAIL);
    return 0;
  }

  long
    del_attachment
    (long given_chain_index, long flag_HEAD_TAIL)
  {
    // it make detachment for subjected chain ends
    // then attach to particles which other chain ends are attached
    mov_attachment(ATTACHED(given_chain_index, (flag_HEAD_TAIL + 1)%2), given_chain_index, flag_HEAD_TAIL);
    return 0;
  }
    
  long
    initial
    (long number_of_chains, long number_of_particles)
  {
    printf("ST:initial:CHAIN_INFORMATION\n");
    N_chains = number_of_chains;
    N_particles = number_of_particles;
    CHAIN = new CHAIN_NODE [N_chains];
    
    for(long i=0; i<N_chains; i++)
      {
        /*
          it allocate the equivalent number of chain ends per particles
        */
        HEAD(i) = i%N_particles;
        TAIL(i) = i%N_particles;
      }
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }

  long
    initial
    (ASSOCIATION& CONNECT)
  {
    N_chains = CONNECT.Nc*CONNECT.Np;
    N_particles = CONNECT.Np;
    CHAIN = new CHAIN_NODE [N_chains];
    long cnt = 0;
    
    for(long i=0; i<N_particles; i++)
      {
        for(long j=i; j<N_particles; j++) // it includes i=j which is loop and prevents j<i means all the pair is once accounted.
          {
            long wij = CONNECT.return_multiplicity(i, j);
            if(i==j)
              {
                // note that weight count chain end rather than chain
                // when i != j, the duplication is removed because we only count the pair with condition j>i
                // when i == j, the duplication happens in the wij.
                wij /=2;
              }
            for(long w=0; w<wij; w++)
              {
                if(cnt>=N_chains)
                  throw std::invalid_argument("connectivity contains too many chains");
                HEAD(cnt) = i;
                TAIL(cnt) = j;
                ++cnt;
              }
          }
      }
    if(cnt!=N_chains)
      throw std::invalid_argument("connectivity chain count is inconsistent");
    printf("end_particle_allocation\n");
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }

  long
    inheritance_chain
    (COND& given_condition, ASSOCIATION& CONNECT)
  {
    N_chains = CONNECT.Nc*CONNECT.Np;
    N_particles = CONNECT.Np;
    // the problem of the original is that 'malloc' will not call the constructor
    CHAIN = new CHAIN_NODE [N_chains];

    long cnt = 0;
    ifstream GIVEN_CHAIN_FILE;
    GIVEN_CHAIN_FILE.open(given_condition("CONTINUATION_CHAIN_FN").c_str());
    string line;
    while(getline(GIVEN_CHAIN_FILE, line))
      cnt ++;
    GIVEN_CHAIN_FILE.clear();

    GIVEN_CHAIN_FILE.seekg(0);
    long N_steps = atoi(given_condition("CONTINUATION_STEP").c_str());
    if(N_steps == -1)
      {
        for(long i=0; i<cnt-1; i++)
          getline(GIVEN_CHAIN_FILE, line);
      }
    else
      {
        for(long i=0; i<N_steps; i++)
          getline(GIVEN_CHAIN_FILE, line);
      }
    
    for(long k=0; k<N_chains; k++)
      GIVEN_CHAIN_FILE >> HEAD(k); // index from 0 to N_chains-1 : HEAD
    
    for(long k=0; k<N_chains; k++)
      GIVEN_CHAIN_FILE >> TAIL(k); // index from N_chains to 2*N_chains-1 : TAIL
    
    GIVEN_CHAIN_FILE.close();
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }
  
  CHAIN_INFORMATION() 
    {
      std::cout << "ERR: CHAIN Class must have initialization argument\n";
    }
  CHAIN_INFORMATION(long number_of_chains, long number_of_particles)
    {
      initial(number_of_chains, number_of_particles);
    }
  CHAIN_INFORMATION(COND& given_condition)
    {
      if(given_condition("tracking_individual_chain") == "TRUE")
        {
          if(given_condition("CONTINUATION_CHAIN") == "TRUE")
            {
              printf("ERR: Without CONTINUATION TOPOLOGICAL INFORMATION (.hash, .weight), the tracking individual chain cannot inherit existing .chain file information\n");
            }
          else
            {
              printf("ST:Constructor:CHAIN_INFORMATION\n");
              initial(atoi(given_condition("N_chains_per_particle").c_str())*atoi(given_condition("Np").c_str()), atoi(given_condition("Np").c_str()));
              printf("END:Constructor:CHAIN_INFORMATION\n");
            }
        }
    }
  CHAIN_INFORMATION(COND& given_condition, ASSOCIATION& CONNECT)
    // it initialize based on the CONNECT information 
    {
      if(given_condition("tracking_individual_chain") == "TRUE")
        {
          if(given_condition("CONTINUATION_CHAIN") == "TRUE")
            {
              inheritance_chain(given_condition, CONNECT);
            }
          else
            {
              printf("ST:Constructor:CHAIN_INFORMATION\n");
              initial(CONNECT);
              printf("END:Constructor:CHAIN_INFORMATION\n");
            }
        }
    }
  virtual ~CHAIN_INFORMATION()
    {
      if(INITIALIZATION)
        delete[] CHAIN;
    }
};



// double inline MICELLE::update_SOFT_REPULSION_P2_particle(POTENTIAL_SET& given_POT, ASSOCIATION& CONNECT, long index_particle)
// {
//   given_POT.force_variables.C_rep_arr(index_particle) = given_POT.force_variables.repulsion_coefficient_base*CONNECT.N_ASSOCIATION_PARTICLE(index_particle)*CONNECT.N_ASSOCIATION_PARTICLE(index_particle);
//   // C_rep(p) = C_0 * p^2
  
//   return given_POT.force_variables.C_rep_arr(index_particle);
// }

#endif
