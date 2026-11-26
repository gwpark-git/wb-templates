#include "../lib/matrix.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>

static void near(double actual, double expected)
{
  if (std::abs(actual - expected) > 1e-12)
    throw std::runtime_error("numeric regression: unexpected result");
}

int main()
{
  MATRIX a(3, 1), b(3, 1);
  a(0) = 1.; a(1) = -2.; a(2) = 3.;
  b(0) = 4.; b(1) = 5.; b(2) = -6.;
  a += b;
  near(a(0), 5.); near(a(1), 3.); near(a(2), -3.);
  matrix_mul(a, -2.);
  near(a(0), -10.); near(a(1), -6.); near(a(2), 6.);
  a.add(b);
  near(a(0), -6.); near(a(1), -1.); near(a(2), 0.);

  // Ownership checks survive optimized builds and exercise repeat reallocation.
  MATRIX owned(2,3,4.); owned.DIAGONALIZATION_INITIAL(); owned.eigen_value[0]=7.;
  MATRIX copied(owned), assigned; assigned=owned; owned=owned;
  near(copied(5),4.); near(assigned.eigen_value[0],7.); near(owned(0),4.);
  copied(0)=9.; copied.eigen_value[0]=3.;
  near(owned(0),4.); near(owned.eigen_value[0],7.);
  for(long i=0; i<500; ++i) {
    owned.initial(3,2,2.); owned.DIAGONALIZATION_INITIAL(); owned.DIAGONALIZATION_INITIAL();
  }
  owned.data_delete(); owned.data_delete();
  if(owned.data || owned.eigen_value || owned.size || owned.INITIALIZATION || owned.DIAGONALIZATION)
    throw std::runtime_error("matrix release must reset ownership");
  bool rejected=false; try { a+=assigned; } catch(const std::invalid_argument&) {rejected=true;}
  if(!rejected) throw std::runtime_error("matrix mismatched += must fail before BLAS");
  rejected=false;try { owned.initial(-1,2); }catch(const std::invalid_argument&) {rejected=true;}
  if(!rejected) throw std::runtime_error("negative matrix shape must be rejected");

  MATRIX output(3,2); output.test_arr();
  std::ostringstream stream; output.fprint(stream);
  std::istringstream parser(stream.str()); double value=0.;
  for(long i=0;i<6;++i) { parser>>value;near(value,i+1.); }
  std::ostringstream skipped;output.fprint_skip(skipped,2);
  std::istringstream skipped_parser(skipped.str());
  for(double expected:{1.,2.,5.,6.}) { skipped_parser>>value;near(value,expected); }
  std::ostringstream transposed;output.fprint_transpose(transposed);
  std::istringstream transposed_parser(transposed.str());
  for(double expected:{1.,3.,5.,2.,4.,6.}) { transposed_parser>>value;near(value,expected); }
  rejected=false;try { output.fprint_skip(stream,0); }catch(const std::invalid_argument&) {rejected=true;}
  if(!rejected) throw std::runtime_error("zero matrix output stride must fail");

  // Strides and untouched guard entries exercise the CBLAS ABI directly.
  double x[] = {1., 99., -2., 99., 3.};
  double y[] = {4., 88., 5., 88., -6.};
  cblas_daxpy(3, 2., x, 2, y, 2);
  near(y[0], 6.); near(y[2], 1.); near(y[4], 0.);
  near(y[1], 88.); near(y[3], 88.);
  near(cblas_ddot(3, x, 2, y, 2), 4.);

  // Association suggestions use a diagonal symmetric-band matrix (k=0).
  const double weights[] = {2., 3., 5.};
  const double probabilities[] = {0.5, 0.25, 0.1};
  double pdf[] = {7., 8., 9.};
  near(cblas_ddot(3, weights, 1, probabilities, 1), 2.25);
  cblas_dsbmv(CblasRowMajor, CblasUpper, 3, 0, 1./2.25,
              weights, 1, probabilities, 1, 0., pdf, 1);
  near(pdf[0], 4./9.); near(pdf[1], 1./3.); near(pdf[2], 2./9.);
  near(pdf[0] + pdf[1] + pdf[2], 1.);

  MATRIX values(3, 1), labels(3, 1);
  values(0) = 3.; values(1) = 1.; values(2) = 2.;
  labels(0) = 30.; labels(1) = 10.; labels(2) = 20.;
  values.sort2(labels);
  near(values(0), 1.); near(values(1), 2.); near(values(2), 3.);
  near(labels(0), 10.); near(labels(1), 20.); near(labels(2), 30.);

  const double start = omp_get_wtime();
  if (omp_get_wtime() < start)
    throw std::runtime_error("non-monotonic wall timer");
  std::cout << "PASS: matrix add/scale, strided axpy, dot, diagonal sbmv, paired sort, timer\n";
  return 0;
}
