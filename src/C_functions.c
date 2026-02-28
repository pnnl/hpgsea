#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


SEXP _C_remove_extreme_gene_sets(const SEXP gene_indices,
                                 const SEXP extreme_set_indices,
                                 const SEXP m) {
  const int n_genes = Rf_length(gene_indices);
  const int n_sets = Rf_length(m);
  const int n_extreme_sets = Rf_length(extreme_set_indices);

  const int *restrict pgene_indices = INTEGER(gene_indices);
  const int *restrict pextreme_set_indices = INTEGER(extreme_set_indices);
  const int *restrict pm = INTEGER(m);

  int length_out = 0;

  for (int j = 0; j < n_extreme_sets; ++j) {
    length_out += pm[pextreme_set_indices[j] - 1];
  }

  length_out = n_genes - length_out;

  SEXP out = PROTECT(Rf_allocVector(INTSXP, length_out));
  int *restrict pout = INTEGER(out);

  int start = 0;
  int end = 0;
  int j = 0;
  int shift_backward = 0;

  for (int i = 0; i < n_sets && j < n_extreme_sets; ++i) {
    start = end;
    end = start + pm[i];

    if (i != (pextreme_set_indices[j] - 1)) {
      for (int k = start; k < end; ++k) {
        pout[k - shift_backward] = pgene_indices[k];
      }
    } else {
      shift_backward += pm[i];
      ++j;
    }
  }

  for (int k = end; k < n_genes; ++k) {
    pout[k - shift_backward] = pgene_indices[k];
  }

  UNPROTECT(1);

  return out;
}


// Calculate enrichment scores for all gene sets
SEXP _C_calc_ES(const SEXP y_prime,
                const SEXP r_prime,
                const SEXP Rsum_ranks,
                const SEXP gene_indices,
                SEXP m,
                const SEXP w,
                const SEXP Rmin_size) {
  const int n_sets = Rf_length(m);

  SEXP ES = PROTECT(Rf_allocVector(REALSXP, n_sets));
  double *restrict pES = REAL(ES);
  memset(pES, 0, n_sets * sizeof(double));

  // This only happens with directional gene sets when the genes from all sets
  // are only up-regulated or only down-regulated.
  if ((int)Rf_length(gene_indices) == 0) {
    UNPROTECT(1);

    return ES;
  }

  const double *restrict py_prime = REAL(y_prime);
  const double *restrict pr_prime = REAL(r_prime);
  const int *restrict pgene_indices = INTEGER(gene_indices);
  int *restrict pm = INTEGER(m);
  const int *restrict pw = INTEGER(w);

  const double sum_ranks = REAL(Rsum_ranks)[0];
  const int min_size = INTEGER(Rmin_size)[0];

  int start = 0;
  int end = 0;

  for (int j = 0; j < n_sets; ++j) {
    // Window of indices for genes in the j-th set
    start = end;
    end = start + pm[j];

    // This only happens with directional gene sets when fewer than min_size
    // genes are up or down.
    if (pm[j] < min_size) {
      pm[j] = 0;

      continue;
    }

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    for (int k = start; k < end; ++k) {
      // Subtract 1 to convert from 1-based to 0-based index
      const int idx = pgene_indices[k] - 1;
      const double r_k = pr_prime[idx];
      const double y_k = py_prime[idx];

      sum_r += r_k;
      sum_y += y_k;
      sum_ry += r_k * y_k;
    }

    // sum_y is 0 if all values for genes in the set are 0
    pES[j] = (sum_y == 0.0) ? NA_REAL :
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
SEXP _C_rep_int(SEXP times, SEXP Rlength_out) {
  const int *restrict ptimes = INTEGER(times);
  const int n_times = Rf_length(times);
  const int length_out = INTEGER(Rlength_out)[0]; // sum(times)

  SEXP out = PROTECT(Rf_allocVector(INTSXP, length_out));
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


// Integer pairing function by Matthew Szudzik. See utils.R for documentation.
SEXP _C_pair_szudzik(const SEXP x,
                     const SEXP y) {
  const int N = Rf_length(x);

  const int *restrict px = INTEGER(x);
  const int *restrict py = INTEGER(y);

  SEXP out = PROTECT(Rf_allocVector(INTSXP, N));
  int *restrict pout = INTEGER(out);

  for (int i = 0; i < N; ++i) {
    pout[i] = px[i] < py[i] ?
      py[i] * py[i] + px[i] :
      px[i] * (px[i] + 1) + py[i];
  }

  UNPROTECT(1);

  return out;
}
