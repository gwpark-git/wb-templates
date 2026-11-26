#ifndef BDMC_HS_MC_H
#define BDMC_HS_MC_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct HS HS;

/* Unit diameter, cubic periodic box, wrapped positions in particle-major order. */
HS *hs_create(int n, double L, uint64_t seed, uint64_t stream, double step,
              int brute, const double *positions);
void hs_destroy(HS *state);
uint64_t hs_sweeps(HS *state, int sweeps);
void hs_copy(const HS *state, double *positions, double *displacements);

/* out: capped minimum squared distance, five cumulative unordered shell
 * counts, and mean squared unwrapped MC displacement. */
void hs_measure(const HS *state, double *out);
int hs_check_lists(const HS *state);
int hs_allowed(const HS *state, int particle, const double *trial);

#ifdef __cplusplus
}
#endif

#endif
