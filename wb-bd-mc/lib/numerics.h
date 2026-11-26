#ifndef NUMERICS_H
#define NUMERICS_H

// GSL supplies the CBLAS routines; OpenMP supplies timing and thread control.
// Simulation indices use standard long, matching the former POSIX index type.
#include <gsl/gsl_cblas.h>
#include <omp.h>

#endif
