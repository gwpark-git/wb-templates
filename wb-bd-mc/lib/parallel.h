
#ifndef PARALLEL_H
#define PARALLEL_H

#define BIT 64
#define TRUE 1
#define FALSE 0

#include <iostream>
#include "handle_association.h"

class LOCK
{
 private:
  long INITIALIZATION;
  
 public:
  bool *locker;
  long NL;
    
  bool&
    operator()
    (long const index_target)
    {
      return locker[index_target];

    }
  bool
    CHECKER
    (long const index_target)
  {
    if(locker[index_target])
      return TRUE;
    return FALSE;

  }
  long
    RESET();
  LOCK();
  LOCK
    (long N_target);
  ~LOCK();
};


#endif
