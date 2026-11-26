

#include "handle_association.h"
#include <vector>
#include <utility>

long
CHAIN_HANDLE::
write
(std::ofstream& file)
{
  long cnt = 0;
  for(long i=0; i<2*N_chains; i++)
    {
      file << CE_ATTACHED_REF(i) << '\t';
    }
  file << std::endl;
  /* FILE1.close(); */
  return 0;
}


long
CHAIN_HANDLE::
mov_attachment
(long target_particle, long given_chain_end_index)
{
  long index_particle = (long)CE_ATTACHED_REF(given_chain_end_index);

  // delete existing information
  long hash_check = get_hash_index(index_particle, given_chain_end_index);
  if(hash_check<0 || target_particle<0 || target_particle>=N_particles)
    throw std::invalid_argument("invalid labelled chain movement");
  for(long i=hash_check; i<P_TOKEN[index_particle]-1; i++)
    {
      PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1];
    }
  // Decrement before clearing: at full capacity the old count is one past
  // the last valid slot. This reverse map must agree with CE_ATTACHED_REF.
  PARTICLE[index_particle][--P_TOKEN[index_particle]] = -1;

  // add new information
  PARTICLE[target_particle][P_TOKEN[target_particle]++] = given_chain_end_index;
  
  // it just change from the existing one to the new one. Therefore, the P_TOKEN will not be changed
  CHAIN_INFORMATION::mov_attachment(target_particle, given_chain_end_index%N_chains, (long)(given_chain_end_index/N_chains)); // it calls existing function in CHAIN_INFORMATION

  return 0;
}


long
CHAIN_HANDLE::
TRACKING_ACTION
(ASSOCIATION& CONNECT, long flag_ACTION, INDEX_MC& IDX)
{
  /*
    It is of importance to using ij indices by identify subjected chain end.
    The basic mechanism in association, handle_association library detaches opponents of subjected chain ends because of convenience. 
    Therefore, the index_subjected_chain_end is using flag_itself and flag_other but the define should be inside sequence of if-phrase.
    It might be changed for future convenience, then it would be changed in the Boltzmann distribution in stochastic_simulation.cpp, ACTION::ADD and ACTION::MOV in handle_association.cpp, del_association and mov_association in association.cpp.
    At this moment, keep the conventional rule even if it is in-convinience with our logic.
    Keep in mind that there is no distingushable for the transition probability of both chain ends belong one chain.
  */
  if(flag_ACTION == INDEX_MC::CANCEL)
    return 0;
  else
    {
      long index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]);
      if(index_subject_chain_end == -1)
        {
          throw std::runtime_error("labelled chains disagree with aggregate topology");
          printf("particle info: (itself, other, new) = (%ld, %ld, %ld)\n", IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_new]);
          printf("connectivity info from itself: ");
          for (long i=0; i<CONNECT.TOKEN[IDX.beads[CONNECT.flag_itself]]; i++)
            printf("(%ld, %ld) = %ld (%ld), ", IDX.beads[CONNECT.flag_itself], i, (long)CONNECT.HASH[IDX.beads[CONNECT.flag_itself]](i), (long)CONNECT.weight[IDX.beads[CONNECT.flag_itself]](i));
          printf("\nconnectivity info from other: ");
          for (long i=0; i<CONNECT.TOKEN[IDX.beads[CONNECT.flag_other]]; i++)
            printf("(%ld, %ld) = %ld (%ld), ", IDX.beads[CONNECT.flag_other], i, (long)CONNECT.HASH[IDX.beads[CONNECT.flag_other]](i), (long)CONNECT.weight[IDX.beads[CONNECT.flag_other]](i));
          
          for(long i=0; i<P_TOKEN[IDX.beads[CONNECT.flag_other]]; i++)
            {
              printf("\n(C, OC)([%ld, %ld] = %ld) = (%ld, %ld)\n", IDX.beads[CONNECT.flag_other], i, PARTICLE[IDX.beads[CONNECT.flag_other]][i], (long)CE_ATTACHED_REF(PARTICLE[IDX.beads[CONNECT.flag_other]][i]), (long)CE_ATTACHED_REF(opp_chain_end_index(PARTICLE[IDX.beads[CONNECT.flag_other]][i])));
            }
        }
	  
      if(flag_ACTION == INDEX_MC::ADD)
        {
          // for loop chains
          // head = subjected chain end, head = other chain end
          mov_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
        }
      else if(flag_ACTION == INDEX_MC::OPP_DEL)
        {
          // for bridge chains
          // head = other chain end, tail = subjected chain end
          // note that it calls the typical 'del' function rather than 'opp_del'
          // the subject of dell is the first argument, which means the first argument will be deleted
	  
          mov_attachment(IDX.beads[CONNECT.flag_itself], index_subject_chain_end);
        }
      else if(flag_ACTION == INDEX_MC::MOV)
        {
          // for bridge chains
          // the following delete scheme is exactly the same for OPP_DEL functionality.
          mov_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
          // Because of the previous detachment (way back to the origin), both of chain ends are attached to the same particle.
          // Therefore, selection HEAD or TAIL is NOT importance to make new bridge.
          // Note that HEAD and TAIL are not related with the degeneracy, which means the checking degeneracy function 'get_index_degeneracy' checked for both pairs for HEAD-TAIL and TAIL-HEAD.
        }
      else
        {
          std::cout << "ERR: No left option in the TRACKING_ACTION function in handle_association.h \n";
        }
    }
  return 0;
}

long
CHAIN_HANDLE::
get_index_degeneracy
(long particle_subject, long particle_target)
{
  /*
    This index function is the core for the tracking algorithm.
    Basically, the tracking individual chain related with the data structure.
    The used adjacence list in our algorithm, it is not distingushable between different chain with the same pair of attachment, but counting multiple connection individually.
    That is one of the reason to reduce computational time dramatically from the original code since we can extracting out the computational overhead for degeneracy.
    To track individual chain, the degeneracy is matter which means we need one more selecting procedure in the degenerated chains.
    This algorithm is selecting one of the degenerated chain using generated random number index.

    It is of importance that the returning chain end index should be related with particle_subjection because the flag controls in the TRACKING_ACTION function is ordered function.
  */
  // H-D07: adjacency stores bridge multiplicity, not chain identity. Select a
  // compatible labeled end uniformly; a loop contributes both distinct ends.
  // Replacing this with one edge per particle pair changes the event measure.
  long degeneracy=0;
  // check the chain: HEAD - particle_subject, TAIL - particle_target
  // note that it include all the chain end attached to particle_subject
  for(long i=0; i<P_TOKEN[particle_subject]; i++)
    {
      if((long)CE_ATTACHED_REF(opp_chain_end_index(PARTICLE[particle_subject][i])) == particle_target)
        // it check the bridges. Therefore, opp_chain_end_index is of importance since PARTICLE[particle_subject][i] is the chain end attached to particle_subject
        degeneracy_index_array[degeneracy++] = PARTICLE[particle_subject][i];
    }
  if (degeneracy < 1)
    {
      printf("ERR: degeneracy is %ld, which is not valid for get_index_degeneracy\n", degeneracy);
      return -1;
    }
  long selecting_index = gsl_rng_uniform_int(r_degeneracy_check, degeneracy);
  long subjected_chain_end = degeneracy_index_array[selecting_index];
  return subjected_chain_end;
}

const long INDEX_MC::N_BOOST_COUNT[] = {0, 2, 2, 3};
/*
  This is defined outside of class definition in the header file.
  The initialization for long CACEL and so on are fine with static constant definition.
  However, in the case of array, it should be defined on here.
  N_BOOST_COUNT[IDX.CANCEL] = 0;
  N_BOOST_COUNT[IDX.ADD] = 2;
  N_BOOST_COUNT[IDX.OPP_DEL] = 2;
  N_BOOST_COUNT[IDX.MOV] = 3;
 */

/*
  the reason for inheritance for reference variable is that the reference variables cannot be along without given variables. In this way, the reference variables automatically allocated with exiting variables with beads array.
*/
INDEX_MC::
INDEX_MC()
{
  initial();
  set_initial_variables();
}

INDEX_MC::
INDEX_MC
(const INDEX_MC& given_IDX)                                                
{
  copy_INDEX(given_IDX);
}

long
INDEX_MC::
initial()
{
  ACTION_ARR[CANCEL] = ACTION::CANCEL ;
  ACTION_ARR[ADD]    = ACTION::ADD    ;
  ACTION_ARR[OPP_DEL]    = ACTION::OPP_DEL;
  ACTION_ARR[MOV]    = ACTION::MOV    ;
  return 0;
}

long
INDEX_MC::
copy_INDEX
(const INDEX_MC& given_IDX)
{
  initial();
  set_initial_variables();
  for(long i=0; i<5; i++)
    beads[i] = given_IDX.beads[i];
  return 0;
}

INDEX_MC&
INDEX_MC::
operator=
(const INDEX_MC& given_IDX)
{
  copy_INDEX(given_IDX);
  return *this;
}

long
INDEX_MC::
set_initial_variables()
{
  // Handling the -1 value as initialization is good part for re-usability from other library.
  for(long i=0; i<5; i++)
    beads[i] = -1; // -1 flag means nothing allocated
  return 0;
}

long
ACTION::
ACT
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, long const IDENTIFIER_ACTION)
{
  IDX.ACTION_ARR[IDENTIFIER_ACTION](POTs, CONNECT, IDX.beads, R_minimum_distance_boost);
  return 0;
}

double
ACTION::
UPDATE_INFORMATION
(ASSOCIATION& CONNECT, INDEX_MC& IDX, long cnt_arr[], long const IDENTIFIER_ACTION)
{
  double time_st = omp_get_wtime();
  for(long i=0; i<IDX.N_BOOST_COUNT[IDENTIFIER_ACTION]; i++)
    {
      CONNECT.update_Z_SUGGESTION_particle(IDX.beads[i]);
      CONNECT.update_dPDF_SUGGESTION_particle(IDX.beads[i]);
      CONNECT.update_dCDF_SUGGESTION_particle(IDX.beads[i]);
    }
  return omp_get_wtime() - time_st;
}

long
ACTION::
CANCEL
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long *index_set, MATRIX* R_minimum_distance_boost)
{
  return 0;
}

long
ACTION::
MOV
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long *index_set, MATRIX* R_minimum_distance_boost)
{
  /*
    It is of importance to nocie that the opp_del_association (i and j) are detaching j chain ends rather than i chain end.
    This is the reason to use 'opp_del' for it's name rather than 'del'.
    In consequence, opp_del_association break the opposite chain ends, then grap this chain end to the subjected particle.
    Some details are described in the opp_del_association_hash file which is the original form in the ASSOCIATION class definition.
    Note that this is related with the Boltzmann distribution generated during stochastic_simulation code is based on the itself particle rather than other particle.
    If this scheme is changed to opposite way, CHAIN_HANDLE class is also affected.
   */
  CONNECT.opp_del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], R_minimum_distance_boost[index_set[CONNECT.flag_itself]](index_set[CONNECT.flag_new]));
  
  return 0;
}

long
ACTION::
OPP_DEL
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long *index_set, MATRIX* R_minimum_distance_boost)
{
  CONNECT.opp_del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  return 0;
}

long
ACTION::
ADD
(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long *index_set, MATRIX* R_minimum_distance_boost)
{
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], R_minimum_distance_boost[index_set[CONNECT.flag_itself]](index_set[CONNECT.flag_new]));
  return 0;
}


long
ACTION::
IDENTIFIER_ACTION_BOOLEAN_BOOST
(ASSOCIATION& CONNECT, INDEX_MC& IDX)
{
  // long &index_itself = index_set[2], &index_other = index_set[0], &index_new = index_set[1];
  // check index_other == index_new
  if (CONNECT.CANCEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      // this is the most probable case for canceling
      return IDX.CANCEL;
    }
  // check index_itself == index_new
  else if (CONNECT.ADD_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_ADD_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.ADD;
      else
        return IDX.CANCEL;
    }
  // check index_itself == index_other
  else if (CONNECT.OPP_DEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_OPP_DEL_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.OPP_DEL;
      else return IDX.CANCEL;
    }
  // with nested case, this is only the case of MOV
  else
    {
      if (CONNECT.CHECK_N_MOV_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.MOV;
      return IDX.CANCEL;
    }
  printf("ERR:TRUTH_MAP::NO SPECIFIED ACTION CASE\n");
  return IDX.CANCEL;
}

long
ACTION::
IDENTIFIER_ACTION_BOOLEAN_BOOST_NO_RESTRICTION
(ASSOCIATION& CONNECT, INDEX_MC& IDX)
{
  // long &index_itself = index_set[2], &index_other = index_set[0], &index_new = index_set[1];
  // check index_other == index_new
  if (TRUTH_MAP::CANCEL_ASSOCIATION_BOOST(CONNECT, IDX.beads))
    {
      // this is the most probable case for canceling
      return IDX.CANCEL;
    }
  // check index_itself == index_new
  else if (TRUTH_MAP::ADD_ASSOCIATION_BOOST(CONNECT, IDX.beads))
    {
      return IDX.ADD;
      // if (TRUTH_MAP::CHECK_N_ADD_ASSOCIATION_BOOST(CONNECT, IDX.beads))
      //   return IDX.ADD;
      // else
      //   return IDX.CANCEL;
    }
  // check index_itself == index_other
  else if (TRUTH_MAP::OPP_DEL_ASSOCIATION_BOOST(CONNECT, IDX.beads))
    {
      return IDX.OPP_DEL;
      // if (TRUTH_MAP::CHECK_N_OPP_DEL_ASSOCIATION_BOOST(CONNECT, IDX.beads))
      //   return IDX.OPP_DEL;
      // else return IDX.CANCEL;
    }
  // with nested case, this is only the case of MOV
  else
    {
      return IDX.MOV;
      // if (TRUTH_MAP::CHECK_N_MOV_ASSOCIATION_BOOST(CONNECT, IDX.beads))
      //   return IDX.MOV;
      // return IDX.CANCEL;
    }
  printf("ERR:TRUTH_MAP::NO SPECIFIED ACTION CASE\n");
  return IDX.CANCEL;
}


long
SEARCHING::
bisection
(MATRIX& given_arr, double p)
{
  long N = given_arr.size;
  long k = N/2;
  // double dk = N/2.;
  long dk = N/2;

  do
    {
      // dk /= 2;
      // the conditional phrase (dk%2 != 0) will return 0 or 1 when the modulo is zero or not, respectively
      // this compensate the loss of approaching because of given domain is in integer rather than real number
      // because of it, the convergence rate is slower than the previous one (log2(Np))
      dk = dk/2 + (dk%2 != 0);
      if (given_arr(k) < p)
        k += dk;
      else if (given_arr(k) > p)
        k -= dk;
      else
        return k;
    } while (dk > 1);
    
  if (given_arr(k) < p)
    k += 1;
  return k;
}

long
SEARCHING::
backtrace
(MATRIX& given_arr, double p)
{
  for(long k= given_arr.size - 1; k >=0; k--)
    {
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return 0;
}

long
SEARCHING::
backtrace_cell_list
(MATRIX& given_arr, long TOKEN, double p, long index_particle, RDIST& R_boost)
{

  for(long k= given_arr.size - 1; k >= given_arr.size - TOKEN; k--)
    {
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return given_arr.size - TOKEN;
}

long
CHAIN_HANDLE::
allocate_array()
{
  PARTICLE = new long* [N_particles];
  P_TOKEN = new long [N_particles];
  
  for(long i=0; i<N_particles; i++)
    {
      /* 
         The PARTICLE is N_particles by N_CE_MAX matrix.
      */
      PARTICLE[i] = new long [N_CE_MAX];
      for(long j=0; j<N_CE_MAX; j++)
        {
          PARTICLE[i][j] = -1; // this means the hash direct no chain ends (index for chain ends is started with 0)
        }
      P_TOKEN[i] = 0;
    }
  return 0;
}

long
CHAIN_HANDLE::
hash_initial
(long seed)
{

  allocate_array();

  for(long i=0; i<N_chains; i++)
    {
      // note that the initial in association.h already allocate HEAD(i) and TAIL(i) as i%N_particles (all connections are loop)
      // hence, this initialization is only allocate the PARTICLE HASH table for further approaches
      PARTICLE[HEAD(i)][P_TOKEN[HEAD(i)]++] = i;
      PARTICLE[TAIL(i)][P_TOKEN[TAIL(i)]++] = i + N_chains;
    }

  // this is for generating random number stream which will be necessary for selecting chain in degenerated state
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r_degeneracy_check = gsl_rng_alloc(T);
  gsl_rng_set(r_degeneracy_check, seed);
  // A loop contributes both endpoints; a particle can own all 2*N_chains ends.
  degeneracy_index_array = new long [2*N_chains];
  return 0;
}

long
CHAIN_HANDLE::
allocate_existing_bridges(ASSOCIATION& CONNECT)
{
  if(!INITIALIZATION || N_particles!=CONNECT.Np || N_chains!=CONNECT.Np*CONNECT.Nc)
    throw std::invalid_argument("chain tracking dimensions do not match connectivity");
  // Validate before mutation. Loop weight counts ends; bridge weight counts chains.
  std::vector<std::pair<long,long>> endpoints;
  endpoints.reserve(N_chains);
  for(long i=0; i<N_particles; ++i) {
    for(long h=0; h<CONNECT.TOKEN[i]; ++h) {
      const long j=static_cast<long>(CONNECT.HASH[i](h));
      const double raw=CONNECT.weight[i](h);
      const long weight=static_cast<long>(raw);
      if(j<0 || j>=N_particles || raw!=weight || weight<0 ||
         (i==j && weight%2) || (i!=j && CONNECT.return_multiplicity(j,i)!=weight))
        throw std::invalid_argument("invalid connectivity in chain reconstruction");
      if(j<i) continue;
      const long count=(i==j ? weight/2 : weight);
      if(count>N_chains-static_cast<long>(endpoints.size()))
        throw std::invalid_argument("too many chains in connectivity");
      for(long w=0; w<count; ++w) endpoints.emplace_back(i,j);
    }
  }
  if(static_cast<long>(endpoints.size())!=N_chains)
    throw std::invalid_argument("chain count in connectivity is inconsistent");
  for(long i=0; i<N_particles; ++i) {
    P_TOKEN[i]=0;
    for(long h=0; h<N_CE_MAX; ++h) PARTICLE[i][h]=-1;
  }
  for(long c=0; c<N_chains; ++c) {
    const long i=endpoints[c].first, j=endpoints[c].second;
    HEAD(c)=i; TAIL(c)=j;
    PARTICLE[i][P_TOKEN[i]++]=c;
    PARTICLE[j][P_TOKEN[j]++]=c+N_chains;
  }
  return 0;
}
