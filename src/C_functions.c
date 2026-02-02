#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>


/*
 * See Rcpp_functions.cpp for documentation.
 */


/* Version 1 ===================================================================
*
* This version is preferred when the CPU cache is large enough to contain every
* true ES. If the cache is too small (which is often the case for older CPUs),
* version 2 will be faster.
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

  double *restrict pES = REAL(ES);
  int *restrict pES_start = INTEGER(ES_start);
  int *restrict pES_end = INTEGER(ES_end);
  int *restrict pES_pos_idx = INTEGER(ES_pos_idx);

  double *restrict py = REAL(y);
  double *restrict pr = REAL(r);
  int *restrict pperm_idx = INTEGER(perm_idx);

  int *restrict pstart_vec = INTEGER(start_vec);
  int *restrict pend_vec = INTEGER(end_vec);
  double *restrict pinv_L2_w = REAL(inv_L2_w);

  // Number of unique gene set sizes
  const int n_sizes = Rf_length(start_vec);

  int offset, start, end, ES_start_i, ES_end_i, ES_pos_idx_i;
  double r_k, y_k, sum_r, sum_y, sum_ry, ES_perm;

  // Loop over permutations
  for (int j = 0; j < batch_size; ++j) {
    offset = j * max_set_size;

    sum_r = 0.0;
    sum_y = 0.0;
    sum_ry = 0.0;

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

      ES_perm = (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_L2_w[i];

      // Update information needed to calculate NES and P-values
      ES_start_i = pES_start[i];
      ES_end_i = pES_end[i];
      ES_pos_idx_i = pES_pos_idx[i];

      if (ES_perm < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm;

        for (int k = ES_start_i; k < ES_pos_idx_i; ++k) {
          pn_as_extreme[k] += ES_perm <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm;

        for (int k = ES_pos_idx_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm >= pES[k];
        }
      }
    } // end unique size loop
  } // end permutation loop

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

  double *restrict pES = REAL(ES);
  int *restrict pES_start = INTEGER(ES_start);
  int *restrict pES_end = INTEGER(ES_end);
  int *restrict pES_pos_idx = INTEGER(ES_pos_idx);

  double *restrict py = REAL(y);
  double *restrict pr = REAL(r);
  int *restrict pperm_idx = INTEGER(perm_idx);

  // Up-regulated
  int *restrict pstart_vec_up = INTEGER(start_vec_up);
  int *restrict pend_vec_up = INTEGER(end_vec_up);
  double *restrict pinv_L3_w_up = REAL(inv_L3_w_up);

  // Down-regulated
  int *restrict pstart_vec_down = INTEGER(start_vec_down);
  int *restrict pend_vec_down = INTEGER(end_vec_down);
  double *restrict pinv_L3_w_down = REAL(inv_L3_w_down);

  // Map from Level 3 to Level 2 to combine directional permutation ES
  int *restrict pmap_L3_to_L2_up = INTEGER(map_L3_to_L2_up);
  int *restrict pmap_L3_to_L2_down = INTEGER(map_L3_to_L2_down);

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

  int offset, start_up, end_up, start_down, end_down;
  int ES_start_i, ES_end_i, ES_pos_idx_i;
  double r_k, y_k, sum_r, sum_y, sum_ry, ES_perm;

  // Loop over permutations
  for (int j = 0; j < batch_size; ++j) {
    offset = j * max_set_size;

    sum_r = 0.0;
    sum_y = 0.0;
    sum_ry = 0.0;

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

      // Reverse loop direction and start at the end because we need to avoid
      // overlap with the values that were selected for the up-regulated genes
      // in the same set.
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
      ES_start_i = pES_start[i];
      ES_end_i = pES_end[i];
      ES_pos_idx_i = pES_pos_idx[i];

      if (ES_perm < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm;

        for (int k = ES_start_i; k < ES_pos_idx_i; ++k) {
          pn_as_extreme[k] += ES_perm <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm;

        for (int k = ES_pos_idx_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm >= pES[k];
        }
      }
    } // end unique up/down pairs loop
  } // end permutation loop

  UNPROTECT(2);

  // Nothing is returned. n_as_extreme, n_perm_neg, sum_perm_pos, and
  // sum_perm_neg are modified in place.
}




// Version 2 ===================================================================

// Not currently available for directional gene sets.

// void C_calc_ES_perm(SEXP n_as_extreme,
//                     SEXP n_perm_neg,
//                     SEXP sum_perm_pos,
//                     SEXP sum_perm_neg,
//                     SEXP ES,
//                     SEXP ES_start,
//                     SEXP ES_end,
//                     SEXP ES_pos_idx,
//                     SEXP y,
//                     SEXP r,
//                     SEXP perm_idx,
//                     const int max_set_size,
//                     const int batch_size,
//                     const double sum_ranks,
//                     SEXP start_vec,
//                     SEXP end_vec,
//                     SEXP inv_L2_w)
// {
//   R_CheckUserInterrupt();
//
//   int *restrict pn_as_extreme = INTEGER(n_as_extreme);
//   int *restrict pn_perm_neg = INTEGER(n_perm_neg);
//   double *restrict psum_perm_pos = REAL(sum_perm_pos);
//   double *restrict psum_perm_neg = REAL(sum_perm_neg);
//
//   double *restrict pES = REAL(ES);
//   int *restrict pES_start = INTEGER(ES_start);
//   int *restrict pES_end = INTEGER(ES_end);
//   int *restrict pES_pos_idx = INTEGER(ES_pos_idx);
//
//   double *restrict py = REAL(y);
//   double *restrict pr = REAL(r);
//   int *restrict pperm_idx = INTEGER(perm_idx);
//
//   int *restrict pstart_vec = INTEGER(start_vec);
//   int *restrict pend_vec = INTEGER(end_vec);
//   double *restrict pinv_L2_w = REAL(inv_L2_w);
//
//   // Number of unique gene set sizes
//   const int n_sizes = Rf_length(start_vec);
//
//   SEXP ES_perm = PROTECT(Rf_allocMatrix(REALSXP, batch_size, n_sizes));
//   double *restrict pES_perm = REAL(ES_perm);
//
//   int offset, start, end;
//   double y_k, r_k, sum_r, sum_y, sum_ry, ES_perm;
//
//   // Loop over permutations
//   for (int j = 0; j < batch_size; ++j) {
//     offset = j * max_set_size;
//
//     sum_y = 0.0;
//     sum_r = 0.0;
//     sum_ry = 0.0;
//
//     // Loop over the unique set sizes
//     for (int i = 0; i < n_sizes; ++i) {
//       start = offset + pstart_vec[i];
//       end = offset + pend_vec[i];
//
//       for (int k = start; k < end; ++k) {
//         y_k = py[pperm_idx[k]];
//         r_k = pr[pperm_idx[k]];
//
//         sum_y += y_k;
//         sum_r += r_k;
//         sum_ry += r_k * y_k;
//       }
//
//       pES_perm[i * batch_size + j] =
//         (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_L2_w[i];
//     } // end unique size loop
//   } // end permutation loop
//
//   // Update information needed to calculate NES and P-values
//   int ES_start_i, ES_end_i, ES_pos_idx_i, n_perm_neg_i;
//   double ES_perm_j, sum_perm_pos_i, sum_perm_neg_i;
//
//   // Loop over unique set sizes
//   for (int i = 0; i < n_sizes; ++i) {
//     ES_start_i = pES_start[i];
//     ES_end_i = pES_end[i];
//     ES_pos_idx_i = pES_pos_idx[i];
//
//     offset = i * batch_size;
//
//     n_perm_neg_i = 0;
//     sum_perm_pos_i = 0.0;
//     sum_perm_neg_i = 0.0;
//
//     // Loop over permutations
//     for (int j = 0; j < batch_size; ++j) {
//       ES_perm_j = pES_perm[offset + j];
//
//       if (ES_perm_j < 0.0) {
//         ++n_perm_neg_i;
//         sum_perm_neg_i -= ES_perm_j;
//
//         for (int k = ES_start_i; k < ES_pos_idx_i; ++k) {
//           pn_as_extreme[k] += ES_perm_j <= pES[k];
//         }
//       } else {
//         sum_perm_pos_i += ES_perm_j;
//
//         for (int k = ES_pos_idx_i; k < ES_end_i; ++k) {
//           pn_as_extreme[k] += ES_perm_j >= pES[k];
//         }
//       }
//     } // end permutation loop
//
//     pn_perm_neg[i] += n_perm_neg_i;
//     psum_perm_neg[i] += sum_perm_neg_i;
//     psum_perm_pos[i] += sum_perm_pos_i;
//   } // end unique size loop
//
//   UNPROTECT(1);
//
//   // Nothing is returned. n_as_extreme, n_perm_neg, sum_perm_pos, and
//   // sum_perm_neg are modified in place.
// }
