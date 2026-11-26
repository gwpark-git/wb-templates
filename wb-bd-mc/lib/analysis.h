#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "matrix.h"

namespace POST_PROCESSING
{
  enum class CorrelationWindow { LegacyFixedWindow, AllOrigins };
  // Existing calls retain the historical fixed window; results are overwritten.
  // AllOrigins uses Nt-lag origins and result.rows lags. center subtracts the
  // full-series mean. Invalid shapes/columns/threads/nonfinite data are rejected.
  double
  compute_autocorrelation_brute_force(MATRIX& time_series_data, const long index,
                                      MATRIX& result, const long index_result,
                                      CorrelationWindow window = CorrelationWindow::LegacyFixedWindow,
                                      bool center = false);

  double
  compute_autocorrelation_brute_force_OMP(MATRIX& time_series_data, const long index,
                              const long N_THREADS,
                              MATRIX& result, const long index_result,
                              CorrelationWindow window = CorrelationWindow::LegacyFixedWindow,
                              bool center = false);

  // Explicit text-table driver; connector mode groups columns by (qx,qy,qz).
  // It does not interpret the undocumented historical .chain stub.
  int autocorrelation_cli(int argc, char* argv[], bool connector_mode);
}

#endif // ANALYSIS_H
