// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>

// This should not be so large that the entire vector of permutation ES does not
// fit in L3 cache
#define BLOCK_SIZE 200


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
                       const float *pES,
                       const int *pES_start,
                       const int *pES_end)
{
  for (int s = 0; s < n_sizes; ++s) {
    int low = pES_start[s];
    const int min_idx = low;
    int high = pES_end[s] - 1;
    pES_pos_idx[s] = pES_end[s];

    while (low <= high) {
      int mid = (low + high) >> 1;

      if (pES[mid] >= 0.0f) {
        if (mid == min_idx || pES[mid - 1] < 0.0f) {
          pES_pos_idx[s] = mid;
          break;
        } else {
          high = mid - 1;
        }
      } else {
        low = mid + 1;
      }
    }
  }
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
//' @param pES_pos_idx pointer to an integer vector; contains the index of the
//'   first positive ES for every unique gene set size.
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
inline void update_output_vectors(int *pn_same_sign,
                                  double *psum_ES_perm,
                                  const int n_sizes,
                                  const int nperm,
                                  const int *pES_start,
                                  const int *pES_end,
                                  const int *pES_pos_idx,
                                  const int *pn_perm_neg,
                                  const float *psum_perm_pos,
                                  const float *psum_perm_neg)
{
  for (int s = 0; s < n_sizes; ++s) {
    const int ES_start_s = pES_start[s];
    const int ES_pos_idx_s = pES_pos_idx[s];
    const int ES_end_s = pES_end[s];
    const int n_perm_neg_s = pn_perm_neg[s];
    const int n_perm_pos_s = nperm - n_perm_neg_s;
    const double sum_perm_neg_s = (double) psum_perm_neg[s];
    const double sum_perm_pos_s = (double) psum_perm_pos[s];

    for (int i = ES_start_s; i < ES_pos_idx_s; ++i) {
      pn_same_sign[i] = n_perm_neg_s;
      psum_ES_perm[i] = sum_perm_neg_s;
    }

    for (int i = ES_pos_idx_s; i < ES_end_s; ++i) {
      pn_same_sign[i] = n_perm_pos_s;
      psum_ES_perm[i] = sum_perm_pos_s;
    }
  }
}


inline void calc_ES_perm_internal(float *pES_perm_vec,
                                  const int n_genes,
                                  const int n_sizes,
                                  const int block_size,
                                  const float *py,
                                  const float *pr,
                                  const int max_size,
                                  const float sum_ranks,
                                  const int *punique_m,
                                  const float *punique_w_float) {
  for (int b = 0; b < block_size; ++b) {
    const Rcpp::IntegerVector rand_idx = dqrng::dqsample_int(
      n_genes,
      max_size
    );
    const int *prand_idx = rand_idx.begin();

    int start = 0;
    int idx = b;
    float sum_r = 0.0f;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int s = 0; s < n_sizes; ++s, idx += block_size) {
      const int end = punique_m[s];

      while (start < end) {
        const int idx_k = prand_idx[start];
        ++start;
        sum_r += pr[idx_k];
        sum_y += py[idx_k];
        sum_ry += pr[idx_k] * py[idx_k];
      }

      pES_perm_vec[idx] = sum_ry / sum_y + (sum_r - sum_ranks) * punique_w_float[s];
    }
  }
}

//' @title Update components for P-values and NES
//'
//' @param pn_as_extreme pointer to an integer vector that stores the number of
//'   permutation ES that are at least as extreme as the ES for each gene set.
//' @param pn_perm_neg pointer to an integer vector that stores the number of
//'   negative permutation ES for each unique gene set size.
//' @param psum_perm_neg pointer to a float vector that stores the absolute sums
//'   of the negative permutation ES.
//' @param psum_perm_pos pointer to a float vector that stores the sums of the
//'   positive permutation ES.
//' @param pES_perm_vec a pointer to a float vector containing permutation ES.
//'   Permutation ES for gene sets of the same size are arranged in blocks of
//'   (at most) 32 consecutive elements.
//' @param pES pointer to a float vector of ES for all gene sets.
//' @param pES_start pointer to an integer vector containing the index of the
//'   first ES for each unique set size.
//' @param pES_pos_idx pointer to an integer vector containing the index of the
//'   first positive ES for each unique set size.
//' @param pES_end pointer to an integer vector containing the 1 more than the
//'   index of the last ES for each unique set size.
//' @param n_sizes integer; the number of unique gene set sizes.
//' @param block_size integer; the size of a block of permutation ES. Between 1
//'   and 32.
//'
//' @details For each unique gene set size, the permutation ES are partitioned
//'   by sign. At the same time, the number of negative permutation ES, as well
//'   as the absolute sums of the negative and positive permutation ES, are
//'   updated. Then, 8 negative ES are compared to a single permutation ES,
//'   iterating over all negative permutation ES before moving to the next 8
//'   negative ES to update n_as_extreme. The process is repeated for the
//'   positive ES and positive permutation ES.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
inline void update_n_as_extreme(int *pn_as_extreme,
                                int *pn_perm_neg,
                                float *psum_perm_neg,
                                float *psum_perm_pos,
                                float *pES_perm_vec,
                                const float *pES,
                                const int *pES_start,
                                const int *pES_pos_idx,
                                const int *pES_end,
                                const int n_sizes,
                                const int block_size) {
  int block_start_perm = 0;

  for (int s = 0; s < n_sizes; ++s) {
    // Partition permutation ES by sign and update the sums ====
    float sum_neg = 0.0f;
    float sum_pos = 0.0f;
    int perm_first_pos = block_start_perm;

    for (int k = 0; k < block_size; ++k) {
      const float swap = pES_perm_vec[block_start_perm + k];

      if (swap < 0.0f) {
        pES_perm_vec[block_start_perm + k] = pES_perm_vec[perm_first_pos];
        pES_perm_vec[perm_first_pos] = swap;
        sum_neg -= swap;
        ++perm_first_pos;
      } else {
        sum_pos += swap;
      }
    }

    perm_first_pos -= block_start_perm; // 0 to block_size
    pn_perm_neg[s] += perm_first_pos;
    psum_perm_neg[s] += sum_neg;
    psum_perm_pos[s] += sum_pos;

    // Comparisons ====
    int i_neg = pES_start[s];
    int i_pos = pES_pos_idx[s];
    int remainder_neg = i_pos - i_neg;
    int remainder_pos = pES_end[s] - i_pos;

    // Negative ----
    while (remainder_neg >= 8) {
      for (int b = 0; b < perm_first_pos; ++b) {
        for (int k = 0; k < 8; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[block_start_perm + b] <= pES[i_neg + k];
        }
      }
      i_neg += 8;
      remainder_neg -= 8;
    }

    if (remainder_neg >= 4) {
      for (int b = 0; b < perm_first_pos; ++b) {
        for (int k = 0; k < 4; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[block_start_perm + b] <= pES[i_neg + k];
        }
      }
      i_neg += 4;
      remainder_neg -= 4;
    }

    if (remainder_neg >= 2) {
      for (int b = 0; b < perm_first_pos; ++b) {
        for (int k = 0; k < 2; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[block_start_perm + b] <= pES[i_neg + k];
        }
      }
      i_neg += 2;
      remainder_neg -= 2;
    }

    if (remainder_neg == 1) {
      for (int b = 0; b < perm_first_pos; ++b) {
        pn_as_extreme[i_neg] += pES_perm_vec[block_start_perm + b] <= pES[i_neg];
      }
    }

    // Positive ----
    while (remainder_pos >= 8) {
      for (int b = perm_first_pos; b < block_size; ++b) {
        for (int k = 0; k < 8; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[block_start_perm + b] >= pES[i_pos + k];
        }
      }
      i_pos += 8;
      remainder_pos -= 8;
    }

    if (remainder_pos >= 4) {
      for (int b = perm_first_pos; b < block_size; ++b) {
        for (int k = 0; k < 4; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[block_start_perm + b] >= pES[i_pos + k];
        }
      }
      i_pos += 4;
      remainder_pos -= 4;
    }

    if (remainder_pos >= 2) {
      for (int b = perm_first_pos; b < block_size; ++b) {
        for (int k = 0; k < 2; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[block_start_perm + b] >= pES[i_pos + k];
        }
      }
      i_pos += 2;
      remainder_pos -= 2;
    }

    if (remainder_pos == 1) {
      for (int b = perm_first_pos; b < block_size; ++b) {
        pn_as_extreme[i_pos] += pES_perm_vec[block_start_perm + b] >= pES[i_pos];
      }
    }

    block_start_perm += block_size;
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
//' @param ES_dbl numeric vector of enrichment scores, sorted in ascending order
//'   by gene set size and then by the values of the ES.
//' @param ES_end integer vector with length equal to \code{n_sizes}. Each
//'   element is 1 more than the index of the last ES for every unique gene set
//'   size.
//' @param y the absolute values of the gene level statistics raised to some
//'   non-negative power \code{alpha}.
//' @param r the ranks of the gene-level statistics.
//' @param max_size integer; the size of the largest gene set.
//' @param Rsum_ranks the sum of the ranks of all gene-level statistics (sum of
//'   the vector \code{r}).
//' @param unique_m integer vector of unique gene set sizes, sorted in ascending
//'   order.
//' @param unique_w integer vector; the differences between the total number of
//'   genes and the elements of \code{unique_m} (number of genes not in the
//'   set).
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
                  const SEXP ES_dbl,
                  const SEXP ES_end,
                  const SEXP y_dbl,
                  const SEXP r_dbl,
                  const int max_size,
                  const SEXP Rsum_ranks,
                  const SEXP unique_m,
                  const SEXP unique_w)
{
  int *pn_same_sign = INTEGER(n_same_sign);
  int *pn_as_extreme = INTEGER(n_as_extreme);
  double *psum_ES_perm = REAL(sum_ES_perm);
  const double *pES_dbl = REAL(ES_dbl);
  const int *pES_end = INTEGER(ES_end);
  const double *py_dbl = REAL(y_dbl);
  const double *pr_dbl = REAL(r_dbl);
  const float sum_ranks = (float) REAL(Rsum_ranks)[0];
  const int *punique_m = INTEGER(unique_m);
  const int *punique_w = INTEGER(unique_w);

  const int n_genes = Rf_length(y_dbl);
  const int n_sets = Rf_length(ES_dbl);
  const int n_sizes = Rf_length(unique_m);

  std::vector<float> ES(n_sets);
  float *pES = &ES[0];
  for (int i = 0; i < n_sets; ++i) {
    pES[i] = (float) pES_dbl[i];
  }

  std::vector<float> y(n_genes);
  std::vector<float> r(n_genes);
  float *py = &y[0];
  float *pr = &r[0];
  for (int i = 0; i < n_genes; ++i) {
    py[i] = (float) py_dbl[i];
    pr[i] = (float) pr_dbl[i];
  }

  std::vector<int> n_perm_neg(n_sizes);
  std::vector<float> sum_perm_neg(n_sizes);
  std::vector<float> sum_perm_pos(n_sizes);
  int *pn_perm_neg = &n_perm_neg[0];
  float *psum_perm_neg = &sum_perm_neg[0];
  float *psum_perm_pos = &sum_perm_pos[0];

  std::vector<int> ES_start(n_sizes);
  int *pES_start = &ES_start[0];
  pES_start[0] = 0;
  for (int i = 1; i < n_sizes; ++i) {
    pES_start[i] = pES_end[i - 1];
  }

  std::vector<int> ES_pos_idx(n_sizes);
  int *pES_pos_idx = &ES_pos_idx[0];
  update_ES_pos_idx(pES_pos_idx, n_sizes, pES, pES_start, pES_end);

  std::vector<float> unique_w_float(n_sizes);
  float *punique_w_float = &unique_w_float[0];
  for (int i = 0; i < n_sizes; ++i) {
    punique_w_float[i] = (float) (1.0 / punique_w[i]);
  }

  std::vector<float> ES_perm_vec(BLOCK_SIZE * n_sizes);
  float *pES_perm_vec = &ES_perm_vec[0];

  dqrng::dqset_seed(seed);
  const int partial_block_size = nperm % BLOCK_SIZE;

  if (partial_block_size != 0) {
    Rcpp::checkUserInterrupt();

    calc_ES_perm_internal(
      pES_perm_vec,
      n_genes,
      n_sizes,
      partial_block_size,
      py,
      pr,
      max_size,
      sum_ranks,
      punique_m,
      punique_w_float
    );

    update_n_as_extreme(
      pn_as_extreme,
      pn_perm_neg,
      psum_perm_neg,
      psum_perm_pos,
      pES_perm_vec,
      pES,
      pES_start,
      pES_pos_idx,
      pES_end,
      n_sizes,
      partial_block_size
    );
  }

  const int n_blocks = nperm / BLOCK_SIZE;

  for (int i = 0; i < n_blocks; ++i) {
    Rcpp::checkUserInterrupt();

    calc_ES_perm_internal(
      pES_perm_vec,
      n_genes,
      n_sizes,
      BLOCK_SIZE,
      py,
      pr,
      max_size,
      sum_ranks,
      punique_m,
      punique_w_float
    );

    update_n_as_extreme(
      pn_as_extreme,
      pn_perm_neg,
      psum_perm_neg,
      psum_perm_pos,
      pES_perm_vec,
      pES,
      pES_start,
      pES_pos_idx,
      pES_end,
      n_sizes,
      BLOCK_SIZE
    );
  }

  update_output_vectors(
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
}


inline void calc_ES_perm_dir_internal(float *pES_perm_vec,
                                      float *pES_perm_down,
                                      float *pES_perm_up,
                                      const int n_genes,
                                      const int n_sizes_down,
                                      const int n_sizes_up,
                                      const int n_pairs,
                                      const int block_size,
                                      const float *py,
                                      const float *pr,
                                      const int max_size,
                                      const float sum_ranks,
                                      const int *punique_m_down,
                                      const int *punique_m_up,
                                      const float *pinv_unique_w_up,
                                      const float *pinv_unique_w_down,
                                      const int *pmap_unique_to_pairs_up,
                                      const int *pmap_unique_to_pairs_down) {
  for (int b = 0; b < block_size; ++b) {
    const Rcpp::IntegerVector rand_idx = dqrng::dqsample_int(
      n_genes,
      max_size
    );
    const int *prand_idx = rand_idx.begin();

    // Up-regulated genes ====
    int start = 0;
    float sum_r = 0.0f;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int i_up = 0; i_up < n_sizes_up; ++i_up) {
      const int end = punique_m_up[i_up];
      const bool empty_set = (start == end);

      while (start < end) {
        const int idx_k = prand_idx[start];
        ++start;
        sum_r += pr[idx_k];
        sum_y += py[idx_k];
        sum_ry += pr[idx_k] * py[idx_k];
      }

      pES_perm_up[i_up] = empty_set ? 0.0f :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_unique_w_up[i_up];
    }

    // Down-regulated genes ====
    start = max_size - 1;
    sum_r = 0.0f;
    sum_y = 0.0f;
    sum_ry = 0.0f;

    for (int i_down = 0; i_down < n_sizes_down; ++i_down) {
      const int end = punique_m_down[i_down];
      const bool empty_set = (start == end);

      // Indices for the down-regulated portions of the gene sets start at the
      // end of rand_idx to avoid selecting the same values that were used for
      // the up-regulated permutation ES.
      while (start > end) {
        const int idx_k = prand_idx[start];
        --start;
        sum_r += pr[idx_k];
        sum_y += py[idx_k];
        sum_ry += pr[idx_k] * py[idx_k];
      }

      pES_perm_down[i_down] = empty_set ? 0.0f :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_unique_w_down[i_down];
    }

    // ES = ES_up - ES_down for each unique pair
    int idx = b;
    for (int s = 0; s < n_pairs; ++s, idx += block_size) {
      pES_perm_vec[idx] =
        pES_perm_up[pmap_unique_to_pairs_up[s]] -
        pES_perm_down[pmap_unique_to_pairs_down[s]];
    }
  }
}


//' @title Permutation Tests for Directional Gene Sets
//'
//' @description Calculate permutation enrichment scores for directional gene
//'   sets and update vectors needed to calculate NES and p-values.
//'
//' @inheritParams calc_ES_perm
//' @param unique_m_up integer vector of the unique number of up-regulated genes
//'   found in the directional gene sets, sorted in ascending order.
//' @param unique_w_up integer vector; the differences between the total number
//'   of genes and the elements of \code{unique_m_up} (unique number of genes
//'   not up-regulated in the set).
//' @param unique_m_down integer vector of the unique number of down-regulated
//'   genes found in the directional gene sets, sorted in ascending order.
//' @param unique_w_down integer vector; the differences between the total
//'   number of genes and the elements of \code{unique_m_down} (unique number of
//'   genes not down-regulated in the set).
//' @param map_unique_to_pairs_up a 1-based integer vector that maps each unique
//'   number of up-regulated genes to the unique pairs of up- and down-regulated
//'   genes.
//' @param map_unique_to_pairs_down a 1-based integer vector that maps each
//'   unique number of down-regulated genes to the unique pairs of up- and
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
                      const SEXP ES_dbl,
                      const SEXP ES_end,
                      const SEXP y_dbl,
                      const SEXP r_dbl,
                      const int max_size,
                      const SEXP Rsum_ranks,
                      const SEXP unique_m_up,
                      const SEXP unique_w_up,
                      const SEXP unique_m_down,
                      const SEXP unique_w_down,
                      SEXP map_unique_to_pairs_up,
                      SEXP map_unique_to_pairs_down)
{
  int *pn_same_sign = INTEGER(n_same_sign);
  int *pn_as_extreme = INTEGER(n_as_extreme);
  double *psum_ES_perm = REAL(sum_ES_perm);
  const double *pES_dbl = REAL(ES_dbl);
  const int *pES_end = INTEGER(ES_end);
  const double *py_dbl = REAL(y_dbl);
  const double *pr_dbl = REAL(r_dbl);
  const float sum_ranks = (float) REAL(Rsum_ranks)[0];
  int *punique_m_up = INTEGER(unique_m_up);
  int *punique_w_up = INTEGER(unique_w_up);
  int *punique_m_down = INTEGER(unique_m_down);
  int *punique_w_down = INTEGER(unique_w_down);
  int *pmap_unique_to_pairs_up = INTEGER(map_unique_to_pairs_up);
  int *pmap_unique_to_pairs_down = INTEGER(map_unique_to_pairs_down);

  const int n_genes = Rf_length(y_dbl);
  const int n_sets = Rf_length(ES_dbl);
  const int n_pairs = Rf_length(map_unique_to_pairs_up); // n_pairs <= n_sets
  const int n_sizes_up = Rf_length(unique_m_up);
  const int n_sizes_down = Rf_length(unique_m_down);

  std::vector<float> ES(n_sets);
  float *pES = &ES[0];
  for (int i = 0; i < n_sets; ++i) {
    pES[i] = (float) pES_dbl[i];
  }

  std::vector<float> y(n_genes);
  std::vector<float> r(n_genes);
  float *py = &y[0];
  float *pr = &r[0];
  for (int i = 0; i < n_genes; ++i) {
    py[i] = (float) py_dbl[i];
    pr[i] = (float) pr_dbl[i];
  }

  // Convert from 1-based to 0-based indices
  for (int i = 0; i < n_pairs; ++i) {
    --pmap_unique_to_pairs_up[i];
    --pmap_unique_to_pairs_down[i];
  }

  std::vector<float> ES_perm_up(n_sizes_up);
  std::vector<float> ES_perm_down(n_sizes_down);
  std::vector<int> n_perm_neg(n_pairs);
  std::vector<float> sum_perm_pos(n_pairs);
  std::vector<float> sum_perm_neg(n_pairs);
  float *pES_perm_up = &ES_perm_up[0];
  float *pES_perm_down = &ES_perm_down[0];
  int *pn_perm_neg = &n_perm_neg[0];
  float *psum_perm_pos = &sum_perm_pos[0];
  float *psum_perm_neg = &sum_perm_neg[0];

  std::vector<int> ES_start(n_pairs);
  int *pES_start = &ES_start[0];
  pES_start[0] = 0;
  for (int i = 1; i < n_pairs; ++i) {
    pES_start[i] = pES_end[i - 1];
  }

  std::vector<int> ES_pos_idx(n_pairs);
  int *pES_pos_idx = &ES_pos_idx[0];
  update_ES_pos_idx(pES_pos_idx, n_pairs, pES, pES_start, pES_end);

  // Update end positions of down-regulated genes to avoid overlap with
  // up-regulated genes
  for (int i = 0; i < n_sizes_down; ++i) {
    punique_m_down[i] = (max_size - 1) - punique_m_down[i];
  }

  std::vector<float> inv_unique_w_up(n_sizes_up);
  float *pinv_unique_w_up = &inv_unique_w_up[0];
  for (int i = 0; i < n_sizes_up; ++i) {
    pinv_unique_w_up[i] = (float) (1.0 / punique_w_up[i]);
  }

  std::vector<float> inv_unique_w_down(n_sizes_down);
  float *pinv_unique_w_down = &inv_unique_w_down[0];
  for (int i = 0; i < n_sizes_down; ++i) {
    pinv_unique_w_down[i] = (float) (1.0 / punique_w_down[i]);
  }

  std::vector<float> ES_perm_vec(BLOCK_SIZE * n_pairs);
  float *pES_perm_vec = &ES_perm_vec[0];

  dqrng::dqset_seed(seed);
  const int partial_block_size = nperm % BLOCK_SIZE;

  if (partial_block_size != 0) {
    Rcpp::checkUserInterrupt();

    calc_ES_perm_dir_internal(
      pES_perm_vec,
      pES_perm_down,
      pES_perm_up,
      n_genes,
      n_sizes_down,
      n_sizes_up,
      n_pairs,
      partial_block_size,
      py,
      pr,
      max_size,
      sum_ranks,
      punique_m_down,
      punique_m_up,
      pinv_unique_w_up,
      pinv_unique_w_down,
      pmap_unique_to_pairs_up,
      pmap_unique_to_pairs_down
    );

    update_n_as_extreme(
      pn_as_extreme,
      pn_perm_neg,
      psum_perm_neg,
      psum_perm_pos,
      pES_perm_vec,
      pES,
      pES_start,
      pES_pos_idx,
      pES_end,
      n_pairs,
      partial_block_size
    );
  }

  const int n_blocks = nperm / BLOCK_SIZE;

  for (int i = 0; i < n_blocks; ++i) {
    Rcpp::checkUserInterrupt();

    calc_ES_perm_dir_internal(
      pES_perm_vec,
      pES_perm_down,
      pES_perm_up,
      n_genes,
      n_sizes_down,
      n_sizes_up,
      n_pairs,
      BLOCK_SIZE,
      py,
      pr,
      max_size,
      sum_ranks,
      punique_m_down,
      punique_m_up,
      pinv_unique_w_up,
      pinv_unique_w_down,
      pmap_unique_to_pairs_up,
      pmap_unique_to_pairs_down
    );

    update_n_as_extreme(
      pn_as_extreme,
      pn_perm_neg,
      psum_perm_neg,
      psum_perm_pos,
      pES_perm_vec,
      pES,
      pES_start,
      pES_pos_idx,
      pES_end,
      n_pairs,
      BLOCK_SIZE
    );
  }

  update_output_vectors(
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
}
