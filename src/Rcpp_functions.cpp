// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>


//' @title Get the Index of the First Positive ES for Every Unique Gene Set Size
//'
//' @description Perform binary searches to determine the index of the first
//'   positive ES for each unique gene set size.
//'
//' @param pES_pos_idx pointer to the integer vector that will hold the results.
//' @param n_sizes integer; the number of unique gene set sizes.
//' @param pES pointer to the numeric vector of enrichment scores, sorted in
//'   ascending order by gene set size and then by the values of the ES.
//' @param pES_start pointer to the integer vector with length equal to
//'   \code{n_sizes}. The vector stores the index of the first ES for every
//'   unique gene set size.
//' @param pES_end pointer to the integer vector with length equal to
//'   \code{n_sizes}. Each element is 1 more than the index of the last ES for
//'   every unique gene set size.
//'
//' @returns Nothing. The vector pointed to by \code{pES_pos_idx} is modified in
//'   place. \code{ES_pos_idx} will contain the index of the first positive ES
//'   for every unique gene set size. Each element is greater than or equal to
//'   the corresponding element of \code{ES_start}. If there are no positive ES
//'   for a particular set size, the index will be the corresponding element of
//'   \code{ES_end}.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
void update_ES_pos_idx(int *pES_pos_idx,
                       const int n_sizes,
                       const double *pES,
                       const int *pES_start,
                       const int *pES_end)
{
  // zero_idx is the minimum index for the ES from gene sets of a particular
  // size. Prevents checking ES from smaller gene sets.
  int zero_idx, low, mid, high;

  for (int i = 0; i < n_sizes; ++i) {
    low = pES_start[i];
    zero_idx = low;
    high = pES_end[i] - 1;
    pES_pos_idx[i] = pES_end[i];

    while (low <= high) {
      mid = (low + high) >> 1;

      if (pES[mid] >= 0.0) {
        if (mid == zero_idx || pES[mid - 1] < 0.0) {
          pES_pos_idx[i] = mid;
          break; // break out of while loop
        } else {
          high = mid - 1;
        }
      } else {
        low = mid + 1;
      }
    }
  }

  // ES_pos_idx is modified in place
}


//' @title Update n_same_sign and sum_ES_perm vectors
//'
//' @description Map the vectors with length equal to the number of unique gene
//'   set sizes to the vectors with length equal to the total number of gene
//'   sets.
//'
//' @param pn_same_sign pointer to an integer vector of zeros with length equal
//'   to the number of enrichment scores. Stores the number of permutation ES
//'   with the same sign as each true ES.
//' @param sum_ES_perm pointer to a numeric vector of zeros with length equal to
//'   the number of enrichment scores. Stores the absolute sum of the
//'   permutation ES with the same sign as each true ES.
//' @param n_sizes integer; number of unique gene set sizes.
//' @param nperm integer; total number of permutations.
//' @param pES_start pointer to an integer vector with length equal to
//'   \code{n_sizes}. Stores the index of the first ES for every unique gene set
//'   size.
//' @param pES_end pointer to an integer vector with length equal to
//'   \code{n_sizes}. Each element is 1 more than the index of the last ES for
//'   every unique gene set size.
//' @param pES_pos_idx pointer to an integer vector; the output of
//'   \code{get_ES_pos_idx()}.
//' @param pn_perm_neg pointer to an integer vector; contains the number of
//'   negative permutation ES for every unique gene set size.
//' @param psum_perm_pos pointer to a numeric vector; contains the sum of the
//'   positive permutation ES for every unique gene set size.
//' @param psum_perm_neg pointer to a numeric vector; contains the absolute sum
//'   of the negative permutation ES for every unique gene set size.
//'
//' @returns Nothing. The vectors \code{n_same_sign} and \code{sum_ES_perm} are
//'   modified in place.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
void update_summary_vectors(int *pn_same_sign,
                            double *psum_ES_perm,
                            const int n_sizes,
                            const int nperm,
                            const int *pES_start,
                            const int *pES_end,
                            const int *pES_pos_idx,
                            const int *pn_perm_neg,
                            const double *psum_perm_pos,
                            const double *psum_perm_neg)
{
  int ES_start_i, ES_end_i, ES_pos_idx_i, n_perm_neg_i, n_perm_pos_i;

  double sum_perm_neg_i, sum_perm_pos_i;

  for (int i = 0; i < n_sizes; ++i) {
    ES_start_i = pES_start[i];
    ES_end_i = pES_end[i];
    ES_pos_idx_i = pES_pos_idx[i];

    n_perm_neg_i = pn_perm_neg[i];
    n_perm_pos_i = nperm - n_perm_neg_i;

    sum_perm_neg_i = psum_perm_neg[i];
    sum_perm_pos_i = psum_perm_pos[i];

    for (int j = ES_start_i; j < ES_pos_idx_i; ++j) {
      pn_same_sign[j] = n_perm_neg_i;
      psum_ES_perm[j] = sum_perm_neg_i;
    }

    for (int j = ES_pos_idx_i; j < ES_end_i; ++j) {
      pn_same_sign[j] = n_perm_pos_i;
      psum_ES_perm[j] = sum_perm_pos_i;
    }
  }
}


//' @title Permutation Tests
//'
//' @description Calculate permutation enrichment scores and update vectors
//'   needed to calculate NES and p-values.
//'
//' @param n_same_sign integer vector of zeros. Used to store the number of
//'   permutation ES with the same sign as each true ES.
//' @param n_as_extreme integer vector of zeros. Used to store the number of
//'   permutation ES that are at least as extreme as each true ES.
//' @param sum_ES_perm numeric vector of zeros. Used to store the absolute sums
//'   of the permutation ES with the same sign as each true ES.
//' @param seed integer or \code{NULL}; seed to obtain reproducible results from
//'   permutation tests.
//' @param nperm integer; total number of permutations.
//' @param ES numeric vector of enrichment scores, sorted in ascending order by
//'   gene set size and then by the values of the ES.
//' @param ES_end integer vector with length equal to \code{n_sizes}. Each
//'   element is 1 more than the index of the last ES for every unique gene set
//'   size.
//' @param y the absolute values of the gene level statistics raised to some
//'   non-negative power \code{alpha}.
//' @param r the ranks of the gene-level statistics.
//' @param max_size integer; the size of the largest gene set.
//' @param sum_ranks the sum of the ranks of all gene-level statistics (sum of
//'   the vector \code{r}).
//' @param L2_m integer vector of unique gene set sizes, sorted in ascending
//'   order.
//' @param L2_w integer vector; the differences between the total number of
//'   genes and the elements of \code{L2_m} (number of genes not in the set).
//'
//' @returns Nothing. The vectors \code{n_same_sign}, \code{n_as_extreme}, and
//'   \code{sum_ES_perm} are modified in place.
//'
//' @seealso
//' https://cran.r-project.org/web/packages/dqrng/vignettes/cpp-api.html
//'
//' @author Tyler Sagendorf
//'
//' @noRd

// [[Rcpp::export(.Cpp_calc_ES_perm)]]
void calc_ES_perm(SEXP n_same_sign,
                  SEXP n_as_extreme,
                  SEXP sum_ES_perm,
                  const Rcpp::Nullable<Rcpp::IntegerVector> seed,
                  const int nperm,
                  const SEXP ES,
                  const SEXP ES_end,
                  const SEXP y,
                  const SEXP r,
                  const int max_size,
                  const double sum_ranks,
                  const SEXP L2_m,
                  const SEXP L2_w)
{
  int *pn_same_sign = INTEGER(n_same_sign);
  int *pn_as_extreme = INTEGER(n_as_extreme);
  double *psum_ES_perm = REAL(sum_ES_perm);

  const double *pES = REAL(ES);
  const int *pES_end = INTEGER(ES_end);

  const double *py = REAL(y);
  const double *pr = REAL(r);
  const int *pend_vec = INTEGER(L2_m);

  // Number of gene-level values
  const int n_genes = Rf_length(y);

  // Number of unique gene set sizes
  const int n_sizes = Rf_length(L2_m);

  int nprotect = 0;

  // Vectors updated with each permutation
  SEXP ES_perm = PROTECT(Rf_allocVector(REALSXP, n_sizes)); ++nprotect;
  double *pES_perm = REAL(ES_perm);

  SEXP n_perm_neg = PROTECT(Rf_allocVector(INTSXP, n_sizes)); ++nprotect;
  int *pn_perm_neg = INTEGER(n_perm_neg);
  memset(pn_perm_neg, 0, n_sizes * sizeof(int));

  SEXP sum_perm_pos = PROTECT(Rf_allocVector(REALSXP, n_sizes)); ++nprotect;
  double *psum_perm_pos = REAL(sum_perm_pos);
  memset(psum_perm_pos, 0, n_sizes * sizeof(double));

  SEXP sum_perm_neg = PROTECT(Rf_allocVector(REALSXP, n_sizes)); ++nprotect;
  double *psum_perm_neg = REAL(sum_perm_neg);
  memset(psum_perm_neg, 0, n_sizes * sizeof(double));

  // Index of the first ES in each unique set size group
  SEXP ES_start = PROTECT(Rf_allocVector(INTSXP, n_sizes)); ++nprotect;
  int *pES_start = INTEGER(ES_start);

  pES_start[0] = 0;

  for (int i = 1; i < n_sizes; ++i) {
    pES_start[i] = pES_end[i - 1];
  }

  // Index of the first positive ES in each unique set size group
  SEXP ES_pos_idx = PROTECT(Rf_allocVector(INTSXP, n_sizes)); ++nprotect;
  int *pES_pos_idx = INTEGER(ES_pos_idx);

  update_ES_pos_idx(pES_pos_idx, n_sizes, pES, pES_start, pES_end);

  // Starting positions for segmented sums
  SEXP start_vec = PROTECT(Rf_allocVector(INTSXP, n_sizes)); ++nprotect;
  int *pstart_vec = INTEGER(start_vec);

  pstart_vec[0] = 0;

  for (int i = 1; i < n_sizes; ++i) {
    pstart_vec[i] = pend_vec[i - 1];
  }

  // Number of genes not in each set. Converted to doubles for permutation ES
  const int *pL2_w = INTEGER(L2_w);
  SEXP L2_w_double = PROTECT(Rf_allocVector(REALSXP, n_sizes)); ++nprotect;
  double *pL2_w_double = REAL(L2_w_double);

  for (int i = 0; i < n_sizes; ++i) {
    pL2_w_double[i] = (double)pL2_w[i];
  }

  dqrng::dqset_seed(seed);

  for (int perm = 0; perm < nperm; ++perm) {
    // Checking for interrupt commands every permutation is extremely slow
    if (perm % 50000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Random sample of max_size integers from 0 to n_genes - 1. Used to select
    // pairs of elements from the vectors y and r.
    const Rcpp::IntegerVector rand_idx = dqrng::dqsample_int(n_genes, max_size);
    const int *prand_idx = rand_idx.begin();

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    for (int i = 0; i < n_sizes; ++i) {
      const int start = pstart_vec[i];
      const int end = pend_vec[i];

      // Segmented sums
      for (int k = start; k < end; ++k) {
        const int idx_k = prand_idx[k];
        const double r_k = pr[idx_k];
        const double y_k = py[idx_k];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      pES_perm[i] = sum_ry / sum_y + (sum_r - sum_ranks) / pL2_w_double[i];
    }

    // Separating the calculation of the ES from updating the summary vectors
    // reduces cache misses.
    for (int i = 0; i < n_sizes; ++i) {
      const double ES_perm_i = pES_perm[i];

      if (ES_perm_i < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm_i;

        const int ES_start_i = pES_start[i];
        const int ES_end_i = pES_pos_idx[i];

        // Iterate over negative ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm_i <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm_i;

        const int ES_start_i = pES_pos_idx[i];
        const int ES_end_i = pES_end[i];

        // Iterate over positive ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm_i >= pES[k];
        }
      }
    }

  }

  // Update n_same_sign and sum_ES_perm
  update_summary_vectors(
    pn_same_sign,
    psum_ES_perm,
    n_sizes,
    nperm,
    pES_start,
    pES_end,
    pES_pos_idx,
    pn_perm_neg,
    psum_perm_pos,
    psum_perm_neg
  );

  UNPROTECT(nprotect);

  // n_same_sign, n_as_extreme, and sum_ES_perm are modified in place
}


//' @title Permutation Tests for Directional Gene Sets
//'
//' @description Calculate permutation enrichment scores for directional gene
//'   sets and update vectors needed to calculate NES and p-values.
//'
//' @inheritParams calc_ES_perm
//' @param L3_m_up integer vector of the unique number of up-regulated genes
//'   found in the directional gene sets, sorted in ascending order.
//' @param L3_w_up integer vector; the differences between the total number of
//'   genes and the elements of \code{L3_m_up} (number of genes not up-regulated
//'   in the set).
//' @param L3_m_down integer vector of the unique number of down-regulated genes
//'   found in the directional gene sets, sorted in ascending order.
//' @param L3_w_down integer vector; the differences between the total number of
//'   genes and the elements of \code{L3_m_down} (number of genes not
//'   down-regulated in the set).
//' @param map_L3_to_L2 a 1-based integer vector that maps each unique number of
//'   up-regulated genes to the unique pairs of up- and down-regulated genes.
//' @param map_L3_to_L2_down a 1-based integer vector that maps each unique
//'   number of down-regulated genes to the unique pairs of up- and
//'   down-regulated genes.
//'
//' @returns Nothing. The vectors \code{n_same_sign}, \code{n_as_extreme}, and
//'   \code{sum_ES_perm} are modified in place.
//'
//' @author Tyler Sagendorf
//'
//' @noRd

// [[Rcpp::export(.Cpp_calc_ES_perm_dir)]]
void calc_ES_perm_dir(SEXP n_same_sign,
                      SEXP n_as_extreme,
                      SEXP sum_ES_perm,
                      const Rcpp::Nullable<Rcpp::IntegerVector> seed,
                      const int nperm,
                      const SEXP ES,
                      const SEXP ES_end,
                      const SEXP y,
                      const SEXP r,
                      const int max_size,
                      const double sum_ranks,
                      const SEXP L3_m_up,
                      const SEXP L3_w_up,
                      const SEXP L3_m_down,
                      const SEXP L3_w_down,
                      SEXP map_L3_to_L2_up,
                      SEXP map_L3_to_L2_down)
{
  // Number of gene-level values
  const int n_genes = Rf_length(y);

  // Number of unique pairs of up- and down-regulated genes
  const int n_pairs = Rf_length(map_L3_to_L2_up);
  const int n_sizes_up = Rf_length(L3_m_up);
  const int n_sizes_down = Rf_length(L3_m_down);

  int *pn_same_sign = INTEGER(n_same_sign);
  int *pn_as_extreme = INTEGER(n_as_extreme);
  double *psum_ES_perm = REAL(sum_ES_perm);

  const double *pES = REAL(ES);
  const int *pES_end = INTEGER(ES_end);

  const double *py = REAL(y);
  const double *pr = REAL(r);

  int *pend_vec_up = INTEGER(L3_m_up);
  int *pend_vec_down = INTEGER(L3_m_down);

  int *pmap_L2_to_L2_up = INTEGER(map_L3_to_L2_up);
  int *pmap_L2_to_L2_down = INTEGER(map_L3_to_L2_down);

  // Convert from 1-based to 0-based indices
  for (int i = 0; i < n_pairs; ++i) {
    --pmap_L2_to_L2_up[i];
    --pmap_L2_to_L2_down[i];
  }

  int nprotect = 0;

  // Vectors updated with each permutation
  SEXP ES_perm_up = PROTECT(Rf_allocVector(REALSXP, n_sizes_up));
  ++nprotect;
  double *pES_perm_up = REAL(ES_perm_up);

  SEXP ES_perm_down = PROTECT(Rf_allocVector(REALSXP, n_sizes_down));
  ++nprotect;
  double *pES_perm_down = REAL(ES_perm_down);

  SEXP n_perm_neg = PROTECT(Rf_allocVector(INTSXP, n_pairs)); ++nprotect;
  int *pn_perm_neg = INTEGER(n_perm_neg);
  memset(pn_perm_neg, 0, n_pairs * sizeof(int));

  SEXP sum_perm_pos = PROTECT(Rf_allocVector(REALSXP, n_pairs)); ++nprotect;
  double *psum_perm_pos = REAL(sum_perm_pos);
  memset(psum_perm_pos, 0, n_pairs * sizeof(double));

  SEXP sum_perm_neg = PROTECT(Rf_allocVector(REALSXP, n_pairs)); ++nprotect;
  double *psum_perm_neg = REAL(sum_perm_neg);
  memset(psum_perm_neg, 0, n_pairs * sizeof(double));

  // Index of the first ES in each unique group of up and down genes
  SEXP ES_start = PROTECT(Rf_allocVector(INTSXP, n_pairs)); ++nprotect;
  int *pES_start = INTEGER(ES_start);

  pES_start[0] = 0;

  for (int i = 1; i < n_pairs; ++i) {
    pES_start[i] = pES_end[i - 1];
  }

  // Index of the first positive ES in each unique set size group
  SEXP ES_pos_idx = PROTECT(Rf_allocVector(INTSXP, n_pairs)); ++nprotect;
  int *pES_pos_idx = INTEGER(ES_pos_idx);

  update_ES_pos_idx(pES_pos_idx, n_pairs, pES, pES_start, pES_end);

  // Starting positions for segmented sums (up-regulated genes)
  SEXP start_vec_up = PROTECT(Rf_allocVector(INTSXP, n_sizes_up));
  ++nprotect;
  int *pstart_vec_up = INTEGER(start_vec_up);

  pstart_vec_up[0] = 0;

  for (int i = 1; i < n_sizes_up; ++i) {
    pstart_vec_up[i] = pend_vec_up[i - 1];
  }

  // Starting positions for segmented sums (down-regulated genes)
  SEXP start_vec_down = PROTECT(Rf_allocVector(INTSXP, n_sizes_down));
  ++nprotect;
  int *pstart_vec_down = INTEGER(start_vec_down);

  pstart_vec_down[0] = 0;

  for (int i = 1; i < n_sizes_down; ++i) {
    pstart_vec_down[i] = pend_vec_down[i - 1];
  }

  // Start from the end to avoid overlap with up-regulated genes
  for (int i = 0; i < n_sizes_down; ++i) {
    pstart_vec_down[i] = (max_size - 1) - pstart_vec_down[i];
    pend_vec_down[i] = (max_size - 1) - pend_vec_down[i];
  }

  // Number of genes not up-regulated and in each set. Converted to doubles for
  // permutation ES
  int *pL3_w_up = INTEGER(L3_w_up);
  SEXP L3_w_up_double = PROTECT(Rf_allocVector(REALSXP, n_sizes_up));
  ++nprotect;
  double *pL3_w_up_double = REAL(L3_w_up_double);

  for (int i = 0; i < n_sizes_up; ++i) {
    pL3_w_up_double[i] = (double)pL3_w_up[i];
  }

  // Number of genes not down-regulated and in each set. Converted to doubles
  // for permutation ES
  int *pL3_w_down = INTEGER(L3_w_down);
  SEXP L3_w_down_double = PROTECT(Rf_allocVector(REALSXP, n_sizes_down));
  ++nprotect;
  double *pL3_w_down_double = REAL(L3_w_down_double);

  for (int i = 0; i < n_sizes_down; ++i) {
    pL3_w_down_double[i] = (double)pL3_w_down[i];
  }

  // Map from the number of unique up- or down-regulated genes to the unique
  // pairs of the number of up- and down-regulated genes. (Level 3 --> Level 2)
  int *pmap_L3_to_L2_up = INTEGER(map_L3_to_L2_up);
  int *pmap_L3_to_L2_down = INTEGER(map_L3_to_L2_down);

  dqrng::dqset_seed(seed);

  for (int perm = 0; perm < nperm; ++perm) {
    // Checking for interrupt commands every permutation is extremely slow
    if (perm % 50000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Random sample of max_size integers from 0 to n_genes - 1. Used to select
    // pairs of elements from the vectors y and r.
    const Rcpp::IntegerVector rand_idx = dqrng::dqsample_int(n_genes, max_size);
    const int *prand_idx = rand_idx.begin();

    double sum_r = 0.0;
    double sum_y = 0.0;
    double sum_ry = 0.0;

    for (int i_up = 0; i_up < n_sizes_up; ++i_up) {
      const int start_up = pstart_vec_up[i_up];
      const int end_up = pend_vec_up[i_up];

      // Segmented sums (up-regulated genes)
      for (int k = start_up; k < end_up; ++k) {
        const int idx_k = prand_idx[k];
        const double r_k = pr[idx_k];
        const double y_k = py[idx_k];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      // If any sets have no up-regulated genes, the first element of
      // pES_perm_up will be 0.
      pES_perm_up[i_up] = (start_up == end_up) ? 0.0 :
        (sum_ry / sum_y) + (sum_r - sum_ranks) / pL3_w_up_double[i_up];
    }

    sum_r = 0.0;
    sum_y = 0.0;
    sum_ry = 0.0;

    for (int i_down = 0; i_down < n_sizes_down; ++i_down) {
      const int start_down = pstart_vec_down[i_down];
      const int end_down = pend_vec_down[i_down];

      // // Segmented sums (down-regulated genes). Indices for the
      // down-regulated portions of the gene sets start at the end of rand_idx
      // to avoid selecting the same values that were used for the up-regulated
      // permutation ES. Note that we only need to avoid overlap for each unique
      // combination of the number of up- and down-regulated genes in a set.
      for (int k = start_down; k > end_down; --k) {
        const int idx_k = prand_idx[k];
        const double r_k = pr[idx_k];
        const double y_k = py[idx_k];

        sum_r += r_k;
        sum_y += y_k;
        sum_ry += r_k * y_k;
      }

      // If any sets have no down-regulated genes, the first element of
      // pES_perm_down will be 0.
      pES_perm_down[i_down] = (start_down == end_down) ? 0.0 :
        (sum_ry / sum_y) + (sum_r - sum_ranks) / pL3_w_down_double[i_down];
    }

    for (int i = 0; i < n_pairs; ++i) {
      // Combine ES_up and ES_down for a single unique pair of the number of up
      // and down-regulated genes. ES = ES_up - ES_down
      const double ES_perm_i =
        pES_perm_up[pmap_L3_to_L2_up[i]] -
        pES_perm_down[pmap_L3_to_L2_down[i]];

      if (ES_perm_i < 0.0) {
        ++pn_perm_neg[i];
        psum_perm_neg[i] -= ES_perm_i;

        const int ES_start_i = pES_start[i];
        const int ES_end_i = pES_pos_idx[i];

        // Iterate over negative ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm_i <= pES[k];
        }
      } else {
        psum_perm_pos[i] += ES_perm_i;

        const int ES_start_i = pES_pos_idx[i];
        const int ES_end_i = pES_end[i];

        // Iterate over positive ES only
        for (int k = ES_start_i; k < ES_end_i; ++k) {
          pn_as_extreme[k] += ES_perm_i >= pES[k];
        }
      }
    }

  }

  // Update n_same_sign and sum_ES_perm
  update_summary_vectors(
    pn_same_sign,
    psum_ES_perm,
    n_pairs,
    nperm,
    pES_start,
    pES_end,
    pES_pos_idx,
    pn_perm_neg,
    psum_perm_pos,
    psum_perm_neg
  );

  UNPROTECT(nprotect);

  // n_same_sign, n_as_extreme, and sum_ES_perm are modified in place
}
