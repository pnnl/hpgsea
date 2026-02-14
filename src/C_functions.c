#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


SEXP _C_rep_int(SEXP times) {
  const int *restrict ptimes = INTEGER(times);

  const int n_times = Rf_length(times);

  int n = 0;

  for (int i = 0; i < n_times; ++i) {
    n += ptimes[i];
  }

  SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
  int *restrict pout = INTEGER(out);

  int val = 0;
  int start = 0;
  int end = 0;

  for (int i = 0; i < n_times; ++i) {
    ++val;
    start = end;
    end = start + ptimes[i];

    for (int j = start; j < end; ++j) {
      pout[j] = val;
    }
  }

  UNPROTECT(1);

  return out;
}


SEXP _C_pairs_not_duplicated(SEXP x, SEXP y) {
  const int N = Rf_length(x);

  if (N == 0) {
    return R_NilValue;
  }

  if (N == 1) {
    return Rf_ScalarInteger(1);
  }

  const int *restrict px = INTEGER(x);
  const int *restrict py = INTEGER(y);

  SEXP temp = PROTECT(Rf_allocVector(INTSXP, N));
  int *restrict ptemp = INTEGER(temp);

  ptemp[0] = 1; // the first pair is not a duplicate, by definition

  int n_unique = 1;

  // Due to how x and y are sorted, duplicate pairs will be contiguous
  for (int i = 1; i < N; ++i) {
    if (!((px[i] == px[i - 1]) & (py[i] == py[i - 1]))) { // not duplicated
      ptemp[n_unique] = i + 1;
      ++n_unique;
    }
  }

  if (n_unique == N) {
    UNPROTECT(1);

    return temp;
  }

  // Subset to the first n_unique elements (indices of non-duplicate pairs)
  SEXP out = PROTECT(Rf_allocVector(INTSXP, n_unique));
  int *restrict pout = INTEGER(out);

  for (int i = 0; i < n_unique; ++i) {
    pout[i] = ptemp[i];
  }

  UNPROTECT(2);

  return out;
}


/*
 * See Rcpp_functions.cpp for documentation.
 */

void C_calc_ES_perm(SEXP n_as_extreme,
                    SEXP n_perm_neg,
                    SEXP sum_perm_pos,
                    SEXP sum_perm_neg,
                    SEXP ES,
                    SEXP ES_start,
                    SEXP ES_end,
                    SEXP ES_pos_idx,
                    SEXP y,
                    SEXP r,
                    SEXP perm_idx,
                    const int max_set_size,
                    const int batch_size,
                    const double sum_ranks,
                    SEXP start_vec,
                    SEXP end_vec,
                    SEXP inv_L2_w)
{
  R_CheckUserInterrupt();

  int *restrict pn_as_extreme = INTEGER(n_as_extreme);
  int *restrict pn_perm_neg = INTEGER(n_perm_neg);
  double *restrict psum_perm_pos = REAL(sum_perm_pos);
  double *restrict psum_perm_neg = REAL(sum_perm_neg);

  const double *restrict pES = REAL(ES);
  const int *restrict pES_start = INTEGER(ES_start);
  const int *restrict pES_end = INTEGER(ES_end);
  const int *restrict pES_pos_idx = INTEGER(ES_pos_idx);

  const double *restrict py = REAL(y);
  const double *restrict pr = REAL(r);
  const int *restrict pperm_idx = INTEGER(perm_idx);

  const int *restrict pstart_vec = INTEGER(start_vec);
  const int *restrict pend_vec = INTEGER(end_vec);
  const double *restrict pinv_L2_w = REAL(inv_L2_w);

  // Number of unique gene set sizes
  const int n_sizes = Rf_length(start_vec);

  int start, end;
  double r_k, y_k;

  SEXP ES_perm_vec = PROTECT(Rf_allocVector(REALSXP, n_sizes));
  double *restrict pES_perm_vec = REAL(ES_perm_vec);

  // Loop over permutations
  for (int j = 0; j < batch_size; ++j) {
    const int offset = j * max_set_size;

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    // Loop over the unique set sizes
    for (int i = 0; i < n_sizes; ++i) {
      start = offset + pstart_vec[i];
      end = offset + pend_vec[i];

      // Update the segmented cumulative sums
      for (int k = start; k < end; ++k) {
        r_k = pr[pperm_idx[k]];
        y_k = py[pperm_idx[k]];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      pES_perm_vec[i] = (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_L2_w[i];
    }

    // It seems silly to store the permutation ES in a vector, only to loop over
    // that vector to retrieve the values shortly after, but separating the
    // calculation of the ES from updating the summary vectors allows for better
    // use of the CPU cache memory.
    for (int i = 0; i < n_sizes; ++i) {
      const double ES_perm = pES_perm_vec[i];

      // Update information needed to calculate NES and P-values
      if (ES_perm < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm;

        const int ES_start_i = pES_start[i];
        const int ES_end_i = pES_pos_idx[i];

        // Iterate over negative ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm;

        const int ES_start_i = pES_pos_idx[i];
        const int ES_end_i = pES_end[i];

        // Iterate over positive ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm >= pES[k];
        }
      }
    } // end unique size loop
  } // end permutation loop

  UNPROTECT(1);

  // Nothing is returned. n_as_extreme, n_perm_neg, sum_perm_pos, and
  // sum_perm_neg are modified in place.
}


// Permutation Tests for Directional Gene Sets ----

void C_calc_ES_perm_dir(SEXP n_as_extreme,
                        SEXP n_perm_neg,
                        SEXP sum_perm_pos,
                        SEXP sum_perm_neg,
                        SEXP ES,
                        SEXP ES_start,
                        SEXP ES_end,
                        SEXP ES_pos_idx,
                        SEXP y,
                        SEXP r,
                        SEXP perm_idx,
                        const int max_set_size,
                        const int batch_size,
                        const double sum_ranks,
                        SEXP start_vec_up,
                        SEXP end_vec_up,
                        SEXP inv_L3_w_up,
                        SEXP start_vec_down,
                        SEXP end_vec_down,
                        SEXP inv_L3_w_down,
                        SEXP map_L3_to_L2_up,
                        SEXP map_L3_to_L2_down)
{
  R_CheckUserInterrupt();

  int *restrict pn_as_extreme = INTEGER(n_as_extreme);
  int *restrict pn_perm_neg = INTEGER(n_perm_neg);
  double *restrict psum_perm_pos = REAL(sum_perm_pos);
  double *restrict psum_perm_neg = REAL(sum_perm_neg);

  const double *restrict pES = REAL(ES);
  const int *restrict pES_start = INTEGER(ES_start);
  const int *restrict pES_end = INTEGER(ES_end);
  const int *restrict pES_pos_idx = INTEGER(ES_pos_idx);

  const double *restrict py = REAL(y);
  const double *restrict pr = REAL(r);
  const int *restrict pperm_idx = INTEGER(perm_idx);

  // Up-regulated
  const int *restrict pstart_vec_up = INTEGER(start_vec_up);
  const int *restrict pend_vec_up = INTEGER(end_vec_up);
  const double *restrict pinv_L3_w_up = REAL(inv_L3_w_up);

  // Down-regulated
  const int *restrict pstart_vec_down = INTEGER(start_vec_down);
  const int *restrict pend_vec_down = INTEGER(end_vec_down);
  const double *restrict pinv_L3_w_down = REAL(inv_L3_w_down);

  // Map from Level 3 to Level 2 to combine directional permutation ES
  const int *restrict pmap_L3_to_L2_up = INTEGER(map_L3_to_L2_up);
  const int *restrict pmap_L3_to_L2_down = INTEGER(map_L3_to_L2_down);

  // Number of unique gene set sizes
  const int n_pairs = Rf_length(map_L3_to_L2_up); // n unique up-down pairs
  const int n_sizes_up = Rf_length(start_vec_up);
  const int n_sizes_down = Rf_length(start_vec_down);

  // Vectors to store permutation ES for the up- and down-regulated genes. The
  // length is equal to the unique number of up- or down-regulated genes,
  // respectively.
  SEXP ES_perm_up = PROTECT(Rf_allocVector(REALSXP, n_sizes_up));
  double *restrict pES_perm_up = REAL(ES_perm_up);

  SEXP ES_perm_down = PROTECT(Rf_allocVector(REALSXP, n_sizes_down));
  double *restrict pES_perm_down = REAL(ES_perm_down);

  int start_up, end_up, start_down, end_down;
  double r_k, y_k, ES_perm;

  // Loop over permutations
  for (int j = 0; j < batch_size; ++j) {
    const int offset = j * max_set_size;

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    // Loop over the unique set sizes (up-regulated)
    for (int i_up = 0; i_up < n_sizes_up; ++i_up) {
      start_up = offset + pstart_vec_up[i_up];
      end_up = offset + pend_vec_up[i_up];

      for (int k = start_up; k < end_up; ++k) {
        r_k = pr[pperm_idx[k]];
        y_k = py[pperm_idx[k]];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      pES_perm_up[i_up] = (start_up == end_up) ? 0.0 :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_L3_w_up[i_up];
    }

    sum_r = 0.0;
    sum_y = 0.0;
    sum_ry = 0.0;

    // Loop over the unique set sizes (down-regulated)
    for (int i_down = 0; i_down < n_sizes_down; ++i_down) {
      start_down = offset + pstart_vec_down[i_down];
      end_down = offset + pend_vec_down[i_down];

      // Indices for the down-regulated portions of the gene sets start at the
      // end of each column of perm_idx to avoid selecting the same values that
      // were used to calculate the ES for up-regulated genes. Note that we only
      // need to avoid overlap for each unique combination of the number of up-
      // and down-regulated genes in a set.
      for (int k = start_down; k > end_down; --k) {
        r_k = pr[pperm_idx[k]];
        y_k = py[pperm_idx[k]];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      pES_perm_down[i_down] = (start_down == end_down) ? 0.0 :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_L3_w_down[i_down];
    }

    for (int i = 0; i < n_pairs; ++i) {
      // Combine ES_up and ES_down for a single unique pair of the number of up
      // and down-regulated genes. ES = ES_up - ES_down
      ES_perm = pES_perm_up[pmap_L3_to_L2_up[i]] -
        pES_perm_down[pmap_L3_to_L2_down[i]];

      // Update information needed to calculate NES and P-values
      if (ES_perm < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm;

        const int ES_start_i = pES_start[i];
        const int ES_end_i = pES_pos_idx[i];

        // Iterate over negative ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm;

        const int ES_start_i = pES_pos_idx[i];
        const int ES_end_i = pES_end[i];

        // Iterate over positive ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm >= pES[k];
        }
      }
    } // end unique up/down pairs loop
  } // end permutation loop

  UNPROTECT(2);

  // Nothing is returned. n_as_extreme, n_perm_neg, sum_perm_pos, and
  // sum_perm_neg are modified in place.
}
