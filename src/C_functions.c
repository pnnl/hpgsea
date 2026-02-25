#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


// Calculate enrichment scores for all gene sets
SEXP _C_calc_ES(const SEXP y_prime,
                const SEXP r_prime,
                const SEXP Rsum_ranks,
                const SEXP i,
                const SEXP m,
                const SEXP w) {
  const int n_sets = Rf_length(m);

  SEXP ES = PROTECT(Rf_allocVector(REALSXP, n_sets));
  double *restrict pES = REAL(ES);
  memset(pES, 0, n_sets * sizeof(double));

  // This only happens with directional gene sets when the genes from all sets
  // are only up-regulated or only down-regulated.
  if ((int)Rf_length(i) == 0) {
    UNPROTECT(1);

    return ES;
  }

  const double *restrict py_prime = REAL(y_prime);
  const double *restrict pr_prime = REAL(r_prime);
  const int *restrict pi = INTEGER(i);
  const int *restrict pm = INTEGER(m);
  const int *restrict pw = INTEGER(w);

  const double sum_ranks = REAL(Rsum_ranks)[0];

  int start = 0;
  int end = 0;

  for (int j = 0; j < n_sets; ++j) {

    // Window of indices for genes in the j-th set
    start = end;
    end = start + pm[j];

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    for (int k = start; k < end; ++k) {
      // Subtract 1 to convert from 1-based to 0-based index
      const int idx = pi[k] - 1;
      const double r_k = pr_prime[idx];
      const double y_k = py_prime[idx];

      sum_r += r_k;
      sum_y += y_k;
      sum_ry += r_k * y_k;
    }

    // Directional gene sets may be empty (all genes in the set are up or down).
    // sum_y is 0 if all values for genes in the set are 0.
    pES[j] = (start == end) ? 0.0 :
      (sum_y == 0.0) ? R_NaReal :
      sum_ry / sum_y + (sum_r - sum_ranks) / (double)pw[j];
  }

  UNPROTECT(1);

  return ES;
}


// Similar to table(), but allows counts to be zero
// group_ids is the output of _C_rep_int()
SEXP _C_group_sizes(const SEXP group_ids,
                    const SEXP n_groups) {
  const int *restrict pgroup_ids = INTEGER(group_ids);

  const int N = Rf_length(group_ids);
  const int NGRP = INTEGER(n_groups)[0];

  SEXP group_sizes = PROTECT(Rf_allocVector(INTSXP, NGRP));
  int *restrict pgroup_sizes = INTEGER(group_sizes);
  memset(pgroup_sizes, 0, NGRP * sizeof(int));

  for (int i = 0; i < N; ++i) {
    const int g = pgroup_ids[i] - 1;
    ++pgroup_sizes[g];
  }

  UNPROTECT(1);

  return group_sizes;
}

// Faster rep.int(seq_along(times), times)
SEXP _C_rep_int(SEXP times) {
  const int *restrict ptimes = INTEGER(times);
  const int n_times = Rf_length(times);

  int n = 0;

  for (int i = 0; i < n_times; ++i) {
    n += ptimes[i];
  }

  SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
  int *restrict pout = INTEGER(out);

  int start = 0;
  int end = 0;

  for (int i = 0; i < n_times; ++i) {
    start = end;
    end = start + ptimes[i];
    const int val = i + 1;

    for (int j = start; j < end; ++j) {
      pout[j] = val;
    }
  }

  UNPROTECT(1);

  return out;
}
