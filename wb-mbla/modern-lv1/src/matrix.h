
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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <limits>
// #include <mkl.h>
#include <math.h>
//#include <omp.h>

extern "C" {
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
}

class MATRIX
{
 public:
  long rows = 0;
  long cols = 0;
  long size = 0;
  double *data = nullptr;  // 1-dimensional array
  bool VECTOR = FALSE; 
  bool INITIALIZATION = FALSE;
  bool DIAGONALIZATION = FALSE;
  double *eigen_value = nullptr;

 private:
  // Vectors own contiguous row-major storage; data is a non-owning GSL view.
  std::vector<double> storage_;
  std::vector<double> eigen_storage_;
 public:
  std::string name;

  long argmax_1d();

  
  long DIAGONALIZATION_INITIAL()                     //@ 
  {
    eigen_storage_.assign(rows, 0.);
    eigen_value = eigen_storage_.data();
    DIAGONALIZATION = TRUE;
    return 0;
  }
  long print_eigen();
  long print_eigen(long n_ele);


  long p_self();
  long initial();
  long initial(long N_r, long N_c);                   
  long initial(long N_r, long N_c, double x);
  long print();
  long print(long n_row, long n_col);

  long fprint_skip(std::ofstream& file, long N_skip);
  long fprint_skip_transpose(std::ofstream& file, long N_skip);
  long fprint_LONG_skip(std::ofstream& file, long N_skip);
  long fprint_LONG_skip_transpose(std::ofstream& file, long N_skip);
  long fprint_LONG_skip_transpose_LIMROWS(std::ofstream& file, long N_skip, long N_lim_rows);
  long fprint(std::ofstream& file);
  long fprint_basic(std::ofstream& file);
  long fprint_transpose(std::ofstream& file);
  long fprint_LONG(std::ofstream& file);
  long fprint_LONG_transpose(std::ofstream& file);
  long fprint_row(std::ofstream& file, long given_row);
  long fprint_out_skip(std::ofstream& file, long N_skip);
  long fprint_out(std::ofstream& file);

  long test_arr();
  long test_arr_symm();
  long nonzero();
  long ABS_cond(double x);
  double ABS_sum();
  long force_clean(double TOLERANCE);
  double sum();
  double average();
  long add(MATRIX& given_MAT);


  long read_exist_data(const char* fn_given_array);

  // Constructor
  MATRIX();
  MATRIX(const char* fn_given_array);
  MATRIX(long N_r, long N_c);
  MATRIX(long N_r, long N_c, double x);
  // Copy Constructor
  long copy_obj(const MATRIX& Mat);
  MATRIX(const MATRIX& Mat);


  // Destructor
  virtual ~MATRIX();

  long data_delete();
  long undefined_error();
  long sort();
  long sort2(MATRIX& given_index);

  inline double norm()
  {
    double result = 0.0;
    for(long  i=0; i<size; i++)
      {
        result += pow(data[i],2.0);
      }
    return sqrt(result);
  }

  
  inline long set_value(double x)
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

  
  inline long index(long i, long j)
  {
    return i*cols+j;
  }
  
  // operator
  // inlined operator overloading
  inline MATRIX& operator=(const MATRIX &Mat)
  {
    copy_obj(Mat);
    return *this;
  }

  inline MATRIX& operator+=(const MATRIX &Mat)
    {
      if (size != Mat.size)
        throw std::invalid_argument("matrix += shape mismatch");
      cblas_daxpy(size, 1.0, Mat.data, 1, data, 1);
      return *this;
    }

  inline double& operator()(long  i, long  j)
  {
    if (i < 0 || j < 0 || i >= rows || j >= cols)
      throw std::out_of_range("matrix row/column index");
    return data[i*cols+j];
  }

  inline double& operator()(long  i)
  {
    if (i < 0 || i >= size) throw std::out_of_range("matrix flat index");
    return data[i];
  }

  
};


inline long make_unit_vector(MATRIX& given_vec)
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

inline long matrix_mul(MATRIX& given_vec, double val)
{
  cblas_dscal(given_vec.size, val, given_vec.data, 1);
  return 0;
}


long nonzero(const MATRIX &A);

long gsl_sol_Ax_b(MATRIX& x_vec, MATRIX& A_mat, MATRIX& b_vec, long ind_st = 0);


#endif 
