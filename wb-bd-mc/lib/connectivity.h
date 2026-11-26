
#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include "numerics.h"
#include "matrix.h"
/* #include "read_file_condition.h" */
#include "file_IO.h"

class CONNECTIVITY
{
 public:
  MATRIX *HASH;
  long *TOKEN;
  
  long Np; // number of particles (particle == bead == micelle)
  long Mc; // maximally avaliable connections
  
  CONNECTIVITY()
    {
      std::cout << "ERR:: Class CONNECTIVITY must have argument\n";
    }

  long
    dynamic_allocation
    (long number_of_particles, long maximum_connections); // for united approach
  CONNECTIVITY
    (long number_of_particles, long maximum_connections);
  CONNECTIVITY
    (COND& given_condition);
  virtual
    ~CONNECTIVITY()
    {
      delete[] HASH;
      delete[] TOKEN;
    }

  long
    read_exist_hash
  (const char* fn_hash, long N_steps);
  
  long
    check_valid
    (long index_particle, long index_target);
};

class CHAIN_NODE
{
  /* 
     This class is designed for suggesting node for connectivity information based on single chain subject.
     Basically, its data structure is given by linked-list, and the index itself is related with the array of object.
  */
 public:
  long HEAD, TAIL; // the design is changed according to the design interface.

  long&
    index
    (long flag_HEAD_TAIL)
    {
      /*
        It prevent the complicate interface design which contains duplicated definition.
      */
      if (flag_HEAD_TAIL == 1)
        return TAIL;
      return HEAD;
      // note that the reference variable can omit the & symbol when return its values.
    }
  CHAIN_NODE()
    { HEAD = -1; TAIL = -1; }
  virtual
    ~CHAIN_NODE(){
  }
};


#endif
