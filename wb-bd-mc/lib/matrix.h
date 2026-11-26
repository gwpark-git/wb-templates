
/*
  #################################################################
  ##   Private MATRIX Library ver 0.3                            ##
  ##   Made by "Gun Woo Park", Email : gunwoo.park@unina.it      ##
  ##   Fellow of SUPOLEN, Marie Curie Actions                    ##
  ##   DICMaPI, University of Naples Federico II                 ##
  #################################################################
*/

// all the integer is converted to long

// @ denoted for tested function
// Lots of functions are based on Call by Reference.
// Notice that for matrix operation, we have to put &.
// Vector can use by rows = 1.

// Note that I seperated library for linear algebra and matrix.
// If some error occurs during compiling, please check to include linear algebra code.

#ifndef MATRIX_H
#define MATRIX_H
#define TRUE 1
#define FALSE 0
#define BIT 64
#include <iostream>
#include <fstream>
#include "numerics.h"
#include <math.h>
#include <stdexcept>
//#include <omp.h>

extern "C" {
#include <gsl/gsl_sort_double.h>
}

class MATRIX
{
 public:
  long rows = 0;
  long cols = 0;
  long size = 0;
  double *data = nullptr;  // owned row-major array
  bool VECTOR = false;
  bool INITIALIZATION = false;
  bool DIAGONALIZATION = false;
  double *eigen_value = nullptr;
  long DIAGONALIZATION_INITIAL()                     //@
  {
    double* next = new double [rows]();
    delete[] eigen_value;
    eigen_value = next;
    DIAGONALIZATION = TRUE;
    return 0;
  }
  long
    print_eigen();
  long
    print_eigen
    (long n_ele);                                //@
  // Public Member Function
  long
    p_self();                                     //@
  long
    initial();
  long
    initial
    (long N_r, long N_c);                    //@
  long
    initial
    (long N_r, long N_c, double x);
  long
    print();                                      //@                     
  long
    print
    (long n_row, long n_col);                  //@

  
  long
    fprint_skip
    (std::ostream& file, long N_skip);
  long
    fprint_skip_transpose
    (std::ostream& file, long N_skip);
  long
    fprint_LONG_skip
    (std::ostream& file, long N_skip);
  long
    fprint_LONG_skip_transpose
    (std::ostream& file, long N_skip);
  long
    fprint_LONG_skip_transpose_LIMROWS
    (std::ostream& file, long N_skip, long N_lim_rows);
  long
    fprint
    (std::ostream& file);                             // ???
  long
    fprint_transpose
    (std::ostream& file);                             // ???
  long
    fprint_LONG
    (std::ostream& file);                             // ???
  long
    fprint_LONG_transpose
    (std::ostream& file);                             // ???
  long
    fprint_row
    (std::ostream& file, long given_row);
  long
    fprint_out_skip
    (std::ostream& file, long N_skip);
  long
    fprint_out
    (std::ostream& file);                             // ???

  long
    test_arr();                                   //@
  long
    test_arr_symm();                              //@
  long
    nonzero();
  long
    ABS_cond
    (double x);
  double
    ABS_sum();
  long
    force_clean(double TOLERANCE);
  double
    sum();
  double
    average();
  long
    add(MATRIX& given_MAT);


  long
  read_exist_data(const char* fn_given_array);

  // Constructor
  MATRIX();                                         //@
  MATRIX
    (long N_r, long N_c);                         //@
  MATRIX
    (long N_r, long N_c, double x);               //@
  // Copy Constructor
  long
    copy_obj
    (const MATRIX& Mat);
  MATRIX
    (const MATRIX& Mat);                        //@
  /* int connect_copy_constructor(const MATRIX& Mat); */
  // Destructor
  virtual
    ~MATRIX();                                        //@

  long
    data_delete();                                //@
  long
    undefined_error();                            //@
  /* long index(long i, long j);                          //@ */ // inlined
  long
    sort();
  long
    sort2
    (MATRIX& given_index);
  /* long sort2(size_t* index); */
  // inlined functions
  // note that function inside class definition is inlined
  // these are frequently called function that have small amount of computation

  double
    norm()
  {
    double result = 0.0;
    for(long  i=0; i<size; i++)
      {
        result += pow(data[i],2.0);
      }
    return sqrt(result);
  }

  
  long
    set_value
    (double x)
  {
    if(INITIALIZATION)
      {
        for(long  i=0; i<size; i++)
          {
            data[i] = x;
          } // i
      }
    else
      {
        undefined_error();
      }
    return 0;
  }

  
  long
    index
    (long i, long j)
  {
    return i*cols+j;
  }
  
  // operator
  // inlined operator overloading
  MATRIX&
    operator=
    (const MATRIX &Mat)
    {
      copy_obj(Mat);
      return *this;
    }

  MATRIX&
    operator+=
    (const MATRIX &Mat)
    {
      // std::cout << "the operator += for MATRIX class is in buggy status\n";
      // std::cout << "if this operator is called, to check the internal problems\n";
      if (rows != Mat.rows || cols != Mat.cols)
        throw std::invalid_argument("MATRIX += requires matching shapes");
      // copy_obj(Mat);
      // for(long i=0; i<size; i++)
      //   {
      //     data[i] += Mat.data[i];
      //   }
      cblas_daxpy(size, 1.0, Mat.data, 1, data, 1);
      return *this;
    }

  double&
    operator()
    (long  i, long  j)
  {
    return data[i*cols+j];
  }

  double&
    operator()
    (long  i)
  {
    return data[i];
  }

  
};


inline long
make_unit_vector
(MATRIX& given_vec)
{
  double norm = 0.;
  for(long  i=0; i<given_vec.size; i++)
    {
      norm += pow(given_vec.data[i], 2.0);
    }
  norm = sqrt(norm);
  for(long  i=0; i<given_vec.size; i++)
    {
      given_vec.data[i] /= norm;
    }
  return norm;
}

inline long
matrix_mul
(MATRIX& given_vec, double val)
{
  cblas_dscal(given_vec.size, val, given_vec.data, 1);
  return 0;
}


long
nonzero
(const MATRIX &A);

#endif 
