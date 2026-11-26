
#include "matrix.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <limits>
#include <memory>

long
MATRIX::
read_exist_data(const char* fn_given_array)
{
  std::ifstream GIVEN_FILE;
  GIVEN_FILE.open(fn_given_array);
  long cnt = 0;
  std::string line;
  long cols = 0;
  while(getline(GIVEN_FILE, line))
    {
      if(cnt ==0)
        {
          std::stringstream is(line);
          double tmp;
          while(is >> tmp)
            {cols ++;}
        }
      cnt ++;
    }
  long rows = cnt;
  GIVEN_FILE.clear();
  GIVEN_FILE.seekg(0);
  
  initial(rows, cols, 0.);
  // printf("rows = %ld, cols = %ld\n", rows, cols);
  
  for(long i=0; i<rows; i++)
    {
      for(long j=0; j<cols; j++)
        {
          GIVEN_FILE >> (*this)(i, j);
          // printf("%ld\t", (*this)(i,j));
        }
      // printf("\n");
    }
  GIVEN_FILE.close();
  return 0;
}

long
nonzero
(const MATRIX &A)
{
  long  cnt=0;
  for(long  i=0; i<A.size; i++)
    {
      if(A.data[i] != 0.0)
        {
          cnt++;
        }
    }
  return cnt;
}



// Constructor
MATRIX::
MATRIX()
{
  initial();  
  //std::cout << "Undefined dimension\n";
}

MATRIX::
MATRIX
(long  N_r, long  N_c)
{
  initial(N_r, N_c);
}


MATRIX::
MATRIX
(long  N_r, long  N_c, double x) // initilization with value x
{
  initial(N_r, N_c);
  set_value(x);
}

// Copy-Constructor
long
MATRIX::
copy_obj
(const MATRIX& Mat)
{
  if (this == &Mat) return 0;
  // H-D08/F11: copied eigendata needs its own allocation. Capture it before
  // initial() releases this object's old storage; trajectory/matrix copy tests
  // mutate each copy independently rather than checking equality alone.
  std::unique_ptr<double[]> values;
  if (Mat.DIAGONALIZATION) {
    values.reset(new double[Mat.rows]);
    for (long i=0; i<Mat.rows; ++i) values[i]=Mat.eigen_value[i];
  }
  initial(Mat.rows, Mat.cols);
  for(long i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
  VECTOR = Mat.VECTOR;
  eigen_value = values.release();
  DIAGONALIZATION = Mat.DIAGONALIZATION;
  return 0;
}

MATRIX::
MATRIX
(const MATRIX& Mat)
{
  copy_obj(Mat);
}

// Destructor
MATRIX::
~MATRIX()
{
  if (INITIALIZATION)
    {
      data_delete();
    }
  if (DIAGONALIZATION)
    {
      delete[] eigen_value;
    }
}


// Public Member Function

long
MATRIX::
print_eigen()
{
  print_eigen(rows);
  return 0;
}

long
MATRIX::
print_eigen
(long  n_ele)
{
  if(DIAGONALIZATION)
    {
      std::cout << "Latent Roots = ";
      for(long  i=0; i<n_ele; i++)
        {
          std::cout << eigen_value[i] << '\t';
        }
      std::cout << std::endl;
    }
  else
    {
      std::cout << "MATRIX DOES NOT CALCULATED LATENT ROOTS\n";
      return -1;
    }
  return 0;
}

long
MATRIX::
print()
{
  print(rows, cols);
  return 0;
}

long
MATRIX::
print
(long  n_row, long  n_col)
{
  long  cnt = 0;
  //std::cout << "MATRIX " << this << std::endl;
  for(long  i=0; i<n_row; i++)
    {
      for(long  j=0; j<n_col; j++)
        {
          printf("%8.5lf\t", data[i*cols+j]);
        }
      printf("\n");
    } // i
  printf("\n");
  return 0;
}

long
MATRIX::
fprint_skip
(std::ostream& file, long  N_skip)
{
  if(N_skip<=0) throw std::invalid_argument("matrix row stride must be positive");
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          file << std::scientific << data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

long
MATRIX::
fprint
(std::ostream& file)
{
  fprint_skip(file, 1);
  return 0;
}

long
MATRIX::
fprint_skip_transpose
(std::ostream& file, long  N_skip)
{
  if(N_skip<=0) throw std::invalid_argument("matrix row stride must be positive");
  long  cnt = 0;
  for(long  j=0; j<cols; j++)
    {
      for(long  i=0; i<rows; i+=N_skip)
        {
          file << std::scientific << data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

long
MATRIX::
fprint_transpose
(std::ostream& file)
{
  fprint_skip_transpose(file, 1);
  return 0;
}

long
MATRIX::
fprint_LONG_skip
(std::ostream& file, long  N_skip)
{
  if(N_skip<=0) throw std::invalid_argument("matrix row stride must be positive");
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          file << (long)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

long
MATRIX::
fprint_LONG
(std::ostream& file)
{
  fprint_LONG_skip(file, 1);
  return 0;
}

long
MATRIX::
fprint_LONG_skip_transpose_LIMROWS
(std::ostream& file, long N_skip, long N_lim_rows)
{
  if(N_lim_rows<0 || N_lim_rows>rows) throw std::out_of_range("matrix output row limit");
  if(N_skip<=0) throw std::invalid_argument("matrix row stride must be positive");
  long  cnt = 0;
  for(long  j=0; j<cols; j++)
    {
      for(long  i=0; i<N_lim_rows; i+=N_skip)
        {
          file << (long)data[index(i, j)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

long
MATRIX::
fprint_LONG_skip_transpose
(std::ostream& file, long  N_skip)
{
  MATRIX::fprint_LONG_skip_transpose_LIMROWS(file, N_skip, rows); // count all rows
  return 0;
}

long
MATRIX::
fprint_LONG_transpose
(std::ostream& file)
{
  fprint_LONG_skip_transpose(file, 1);
  return 0;
}


long
MATRIX::
fprint_row
(std::ostream& file, long  given_row)
{
  if(given_row<0 || given_row>=rows) throw std::out_of_range("matrix row index");
  long  cnt = 0;
  long  i=given_row;
  for(long  j=0; j<cols; j++)
    {
      file << std::scientific << data[index(i, j)] << '\t';
    }
  file << std::endl;
  return 0;
}

long
MATRIX::
fprint_out_skip
(std::ostream& file, long  N_skip)
{
  if(N_skip<=0) throw std::invalid_argument("matrix row stride must be positive");
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          file << std::scientific << data[index(i,j)] << '\t';
        }
      file << std::endl;
    } // i
  file << std::endl;
  return 0;

}

long
MATRIX::
fprint_out
(std::ostream& file)
{
  fprint_out_skip(file, 1);
  return 0;
}

long
MATRIX::
p_self()
{
  std::cout << this;
  return 0;
}

long
MATRIX::
nonzero()
{
  for(long  i=0; i<size; i++)
    {
      if(data[i]!=0)
        {
          printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
        }
    }
  return 0;
}

long
MATRIX::
ABS_cond
(double x)
{
  for(long  i=0; i<size; i++)
    {
      if(fabs(data[i]) > fabs(x))
        {
          printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",(int)i/cols,(int)i%cols, data[i]);
        }
    }
  return 0;
}

// Private Member Function
long
MATRIX::
initial()
{
  delete[] data;
  delete[] eigen_value;
  INITIALIZATION = FALSE;
  DIAGONALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  eigen_value = NULL;
  VECTOR = FALSE;
  return 0;
}

long
MATRIX::
initial
(long  N_r, long  N_c)
{
  if (N_r < 0 || N_c < 0 || (N_c && N_r > std::numeric_limits<long>::max()/N_c))
    throw std::invalid_argument("invalid MATRIX shape");
  double* next = new double[N_r*N_c]();
  delete[] data;
  delete[] eigen_value;
  eigen_value = nullptr;
  INITIALIZATION = TRUE;
  DIAGONALIZATION = FALSE;
  rows = N_r;
  cols = N_c;
  size = rows*cols;
  VECTOR = (N_r == 1 || N_c == 1);
  data = next;
  return 0;
}

long
MATRIX::
initial
(long  N_r, long  N_c, double x)
{
  initial(N_r, N_c);
  for(long  i=0; i<size; i++)
    {
      data[i] = x;
    }
  return 0;
}


long
MATRIX::
data_delete()
{
  return initial();
}

long
MATRIX::
undefined_error()
{
  std::cout << "MATRIX is NOT constructed\n" << std::endl;
  return -1;
}


long
MATRIX::
test_arr()
{
  for(long  i=0; i<size; i++)
    {
      data[i] = i+1;
    } // i
  return 0;
}

long
MATRIX::
test_arr_symm()
{
  for(long  i=0; i<rows; i++)
    {
      for(long  j=0; j<cols; j++)
        {
          data[i*cols+j] = 1.0/(1.0+i+j);
        }
    }
  return 0;
}

double
MATRIX::
ABS_sum()
{
  double result = 0.0;
  double tmp_x;
  for(long  i=0; i<size; i++)
    {
      tmp_x = data[i];
      if(tmp_x < 0)
        result -= tmp_x;
      else
        result += tmp_x;
    }
  return result;
}

double
MATRIX::
sum()
{
  double result = 0.0;
  for(long  i=0; i<size; i++)
    {
      result += data[i];
    }
  return result;
}

double
MATRIX::
average()
{
  double result = sum();
  return result/(double)size;
}

long
MATRIX::
add
(MATRIX& given_MAT)
{
  if (size != given_MAT.size)
    {
      std::cout << "ERR: MATRIX::add different size\n";
      return -1;
    }
  cblas_daxpy(size, 1.0, given_MAT.data, 1, data, 1);
  return 0;
}

long
MATRIX::
sort()
{
  gsl_sort(data, 1, size);
  return 0;
}


// long MATRIX::sort2(MATRIX& index)
// {

// long MATRIX::sort2(MATRIX& index)
// {
//   gsl_sort2(data, 1, index.data, 1, size);
//   return 0;
// }

// long MATRIX::sort2(size_t* index)
// {
//   // it is of importance that the gsl_sort2 is not supported by the SCOPE_GRID server.
//   // hence, we are applying new scheme
//   // gsl_sort2(data, 1, index.data, 1, size);
  
//   gsl_sort_index(index, data, 1, size);
//   gsl_sort(data, 1, size);
//   return 0;
// }

long
MATRIX::
sort2
(MATRIX& index)
{

  gsl_sort2(data, 1, index.data, 1, size);
  return 0;
}
