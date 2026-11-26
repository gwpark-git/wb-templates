
#include "matrix.h"


long
MATRIX::
argmax_1d()
{
  if (size == 0) throw std::invalid_argument("argmax of empty matrix");
  long argmax=0;
  double max = data[argmax];
  for (long i=1; i<size; i++)
    {
      if (data[i] > max)
	{
	  max = data[i];
	  argmax = i;
	}
    }
  return argmax;
}

long
MATRIX::
read_exist_data(const char* fn_given_array)
{
  // Reject ragged, nonfinite and partly parsed tables before exposing a GSL view.
  std::ifstream input(fn_given_array);
  if (!input) throw std::runtime_error(std::string("cannot open table: ") + fn_given_array);
  std::vector<double> values;
  std::string line;
  long nrows = 0, ncols = 0;
  while (std::getline(input, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.resize(comment);
    std::istringstream row(line);
    std::vector<double> fields;
    double v;
    while (row >> v) {
      if (!std::isfinite(v)) throw std::runtime_error("nonfinite table value");
      fields.push_back(v);
    }
    if (!row.eof()) throw std::runtime_error("invalid numeric table token");
    if (fields.empty()) continue;
    if (ncols == 0) ncols = fields.size();
    if (static_cast<long>(fields.size()) != ncols) throw std::runtime_error("ragged table");
    ++nrows;
    values.insert(values.end(), fields.begin(), fields.end());
  }
  if (nrows == 0) throw std::runtime_error("empty numeric table");
  initial(nrows, ncols);
  std::copy(values.begin(), values.end(), data);
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
}

MATRIX::
MATRIX
(const char* fn_given_array)
{
  read_exist_data(fn_given_array);
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
  if (this == &Mat) return 0; // Self-assignment must preserve owned values.
  initial(Mat.rows, Mat.cols);
  for(long i=0; i<size; i++)
    {
      data[i] = Mat.data[i];
    }
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
~MATRIX() = default;


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
(std::ofstream& file, long  N_skip)
{
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          file << std::scientific << data[index(j, i)] << '\t';
        }
      file << std::endl;
    } // i
  return 0;
}

long
MATRIX::
fprint_basic
(std::ofstream& file)
{
  for(long i=0; i<rows; i++)
    {
      for(long j=0; j<cols; j++)
	{
	  file << std::scientific << data[index(i,j)] << '\t';
	} //j
      file << std::endl;
    } //i
  return 0;
}


long
MATRIX::
fprint
(std::ofstream& file)
{
  fprint_skip(file, 1);
  return 0;
}

long
MATRIX::
fprint_skip_transpose
(std::ofstream& file, long  N_skip)
{
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
(std::ofstream& file)
{
  fprint_skip_transpose(file, 1);
  return 0;
}

long
MATRIX::
fprint_LONG_skip
(std::ofstream& file, long  N_skip)
{
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
(std::ofstream& file)
{
  fprint_LONG_skip(file, 1);
  return 0;
}

long
MATRIX::
fprint_LONG_skip_transpose_LIMROWS
(std::ofstream& file, long N_skip, long N_lim_rows)
{
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
(std::ofstream& file, long  N_skip)
{
  MATRIX::fprint_LONG_skip_transpose_LIMROWS(file, N_skip, rows); // count all rows
  return 0;
}

long
MATRIX::
fprint_LONG_transpose
(std::ofstream& file)
{
  fprint_LONG_skip_transpose(file, 1);
  return 0;
}


long
MATRIX::
fprint_row
(std::ofstream& file, long  given_row)
{
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
(std::ofstream& file, long  N_skip)
{
  long  cnt = 0;
  for(long  i=0; i<rows; i+=N_skip)
    {
      for(long  j=0; j<cols; j++)
        {
          file << std::scientific << data[cnt++] << '\t';
        }
      file << std::endl;
    } // i
  file << std::endl;
  return 0;

}

long
MATRIX::
fprint_out
(std::ofstream& file)
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
          printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",i/cols,i%cols, data[i]);
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
          printf("nonzero (i,j) = (%ld, %ld) = %6.3e\n",i/cols,i%cols, data[i]);
        }
    }
  return 0;
}

// Private Member Function
long
MATRIX::
initial()
{
  INITIALIZATION = FALSE;
  DIAGONALIZATION = FALSE;
  rows = 0;
  cols = 0;
  size = 0;
  data = NULL;
  name = "NOT_DEFINED";
  return 0;
}

long
MATRIX::
initial
(long  N_r, long  N_c)
{

  if (N_r < 0 || N_c < 0 || (N_c && N_r > std::numeric_limits<long>::max()/N_c))
    throw std::invalid_argument("invalid matrix dimensions");
  rows = N_r;
  cols = N_c;
  size = rows*cols;
  storage_.assign(size, 0.);
  data = storage_.data();
  eigen_storage_.clear();
  eigen_value = nullptr;
  INITIALIZATION = TRUE;
  DIAGONALIZATION = FALSE;
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
  storage_.clear();
  INITIALIZATION = FALSE;
  data = nullptr;
  rows = cols = size = 0;
  return 0;
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



long
MATRIX::
sort2
(MATRIX& index)
{

  gsl_sort2(data, 1, index.data, 1, size);
  return 0;
}


// Solve through GSL's row-major view. LU overwrites the input matrix storage.
long
gsl_sol_Ax_b(MATRIX& x_vec, MATRIX& A_mat, MATRIX& b_vec, long ind_st) // ind_st is default value of 0
{
    // Note that the self-developed MATRIX class has nice compatibility with gsl_vector or gsl_matrix.
    // Therefore, we can simply claim for the mastrix and vector like below codes.
    // Note that those view functions does not make a new array.
    // Which means A_gsl and b_gsl are not necessary to be free.

  long N = b_vec.size;
  if (N == 0 || A_mat.rows != N || A_mat.cols != N || ind_st < 0 || x_vec.size < N + ind_st)
    throw std::invalid_argument("LU solve dimensions");
  gsl_matrix_view A_gsl = gsl_matrix_view_array(A_mat.data, N, N);
  gsl_vector_view b_gsl = gsl_vector_view_array(b_vec.data, N);
  gsl_vector *x = gsl_vector_alloc(N);

  int s;
  gsl_permutation *p = gsl_permutation_alloc(N);
  // GSL overwrites A with LU factors. The pressure assembly rewrites its band
  // before every solve; do not reuse this matrix as an untouched coefficient table.
  int status = gsl_linalg_LU_decomp(&A_gsl.matrix, p, &s);
  if (status == GSL_SUCCESS) status = gsl_linalg_LU_solve(&A_gsl.matrix, p, &b_gsl.vector, x);
  if (status != GSL_SUCCESS) {
    gsl_permutation_free(p);
    gsl_vector_free(x);
    throw std::runtime_error("GSL LU solve failed");
  }

  for (long i=0; i<N; i++)
    {
      x_vec(i + ind_st) = gsl_vector_get(x, i);
    }
  gsl_permutation_free(p);
  gsl_vector_free(x);
  return 0;
}
