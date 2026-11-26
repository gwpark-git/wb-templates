
#include "connectivity.h"

long
CONNECTIVITY::
dynamic_allocation
(long number_of_particles, long maximum_connections)
{
  Np = number_of_particles;
  Mc = maximum_connections;

  HASH = new MATRIX [Np];
  TOKEN = new long [Np];
  for(long i=0; i<Np; i++)
    {
      HASH[i].initial(Mc, 1, -1);
      HASH[i](0) = i; // it will be ignored when it inherit from the given index table
      TOKEN[i] = 1;
    }
  return 0;
}

CONNECTIVITY::
CONNECTIVITY
(COND& given_condition)
{
  dynamic_allocation(atol(given_condition("Np").c_str()), 5 + 2*atol(given_condition("N_chains_per_particle").c_str()) + atol(given_condition("tolerance_allowing_connections").c_str())); // 5 is added for preventing overwhelemd cases
  if (given_condition("CONTINUATION_CONNECTION")=="TRUE")
    {
      read_exist_hash(given_condition("CONTINUATION_HASH_FN").c_str(), atoi(given_condition("CONTINUATION_STEP").c_str()));
    }
  else
    {
      for(long i=0; i<Np; i++)
        {
          HASH[i](0) = i;
        }
    }
}


CONNECTIVITY::
CONNECTIVITY
(long number_of_particles, long maximum_connections)
{
  dynamic_allocation(number_of_particles, maximum_connections);
}

long
CONNECTIVITY::
read_exist_hash
(const char* fn_hash, long N_steps)
{
  ifstream GIVEN_HASH;
  GIVEN_HASH.open(fn_hash);
  long cnt = 0;
  string line;
  while(getline(GIVEN_HASH, line))
    {
      cnt++;
    }
  GIVEN_HASH.clear();
  GIVEN_HASH.seekg(0);
  if(N_steps == -1)
    {
      for(long i=0; i<cnt-Np; i++)
        {
          getline(GIVEN_HASH, line); 
        }
    }
  else
    {
      for(long i=0; i<N_steps*Np; i++)
        {
          getline(GIVEN_HASH, line);
        }
    }
  // flag will set for the last time steps
  // From here, the Np lines will left as initial conditions
  string str;
  long i=0;
  long count_lines = 0;
  while(count_lines++ < Np)
        // getline(GIVEN_HASH, str))
    {
      getline(GIVEN_HASH, str);
      istringstream ss(str);
      long hash_index;
      long tmp_cnt = 0;
      while(ss >> hash_index)
        {
          HASH[i](tmp_cnt) = hash_index;
          tmp_cnt++;
        }
      TOKEN[i] = tmp_cnt;
      i++;
    }
  GIVEN_HASH.close();
  return 0;
}

long
CONNECTIVITY::
check_valid
(long index_particle, long index_target)
{
  if (index_target > TOKEN[index_particle])
    {
      std::cout << "ERR:: Class CONNECTIVITY:: index_target must larger then TOKEN[index_particle]\n";
    }
  return TRUE;
}

