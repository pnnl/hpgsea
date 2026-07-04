// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>

#define BLOCK_SIZE 200


// Gene-level values
typedef struct {
  float r;
  float y;
} gene_data_t;

// For a given set size, defines the index of the first ES, the index of the
// first positive ES, and the index of the last ES + 1.
typedef struct {
  int start;
  int first_positive;
  int end;
} ES_bounds_t;

// For a given set size, count the number of negative permutation ES and compute
// the absolute sums of the negative and positive ES, separately.
typedef struct {
  int n_negative;
  float sum_negative;
  float sum_positive;
} perm_stats_t;

typedef struct {
  int up;
  int down;
} pair_map_t;


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
void update_first_positive(ES_bounds_t *pES_bounds,
                           const int n_sizes,
                           const float *pES)
{
  for (int s = 0; s < n_sizes; ++s) {
    int low = pES_bounds[s].start;
    const int min_idx = low;
    int high = pES_bounds[s].end - 1;
    pES_bounds[s].first_positive = pES_bounds[s].end;

    while (low <= high) {
      int mid = (low + high) >> 1;

      if (pES[mid] >= 0.0f) {
        if (mid == min_idx || pES[mid - 1] < 0.0f) {
          pES_bounds[s].first_positive = mid;
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
                                  const ES_bounds_t *pES_bounds,
                                  const perm_stats_t *pperm_stats)
{
  for (int s = 0; s < n_sizes; ++s) {
    const int start = pES_bounds[s].start;
    const int first_positive = pES_bounds[s].first_positive;
    const int end = pES_bounds[s].end;

    const int n_negative = pperm_stats[s].n_negative;
    const int n_positive = nperm - n_negative;
    const double sum_negative = (double) pperm_stats[s].sum_negative;
    const double sum_positive = (double) pperm_stats[s].sum_positive;

    for (int i = start; i < first_positive; ++i) {
      pn_same_sign[i] = n_negative;
      psum_ES_perm[i] = sum_negative;
    }

    for (int i = first_positive; i < end; ++i) {
      pn_same_sign[i] = n_positive;
      psum_ES_perm[i] = sum_positive;
    }
  }
}


inline void calc_ES_perm_internal(float *pES_perm_vec,
                                  const int n_genes,
                                  const int n_sizes,
                                  const int block_size,
                                  const gene_data_t *pgene_data,
                                  const int max_size,
                                  const float sum_ranks,
                                  const int *punique_m,
                                  const float *pinv_w) {
  for (int b = 0; b < block_size; ++b) {
    const Rcpp::IntegerVector random_indices = dqrng::dqsample_int(
      n_genes,
      max_size
    );
    const int *prandom_indices = random_indices.begin();

    int perm_idx = b;
    int start = 0;
    float sum_r = 0.0f;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int s = 0; s < n_sizes; ++s, perm_idx += block_size) {
      const int end = punique_m[s];

      while (start < end) {
        const int rand_idx = prandom_indices[start++];
        sum_r += pgene_data[rand_idx].r;
        sum_y += pgene_data[rand_idx].y;
        sum_ry += pgene_data[rand_idx].r * pgene_data[rand_idx].y;
      }

      pES_perm_vec[perm_idx] = (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_w[s];
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
//'   (at most) BLOCK_SIZE consecutive elements.
//' @param pES pointer to a float vector of ES for all gene sets.
//' @param pES_start pointer to an integer vector containing the index of the
//'   first ES for each unique set size.
//' @param pES_pos_idx pointer to an integer vector containing the index of the
//'   first positive ES for each unique set size.
//' @param pES_end pointer to an integer vector containing the 1 more than the
//'   index of the last ES for each unique set size.
//' @param n_sizes integer; the number of unique gene set sizes.
//' @param block_size integer; the size of a block of permutation ES. Between 1
//'   and BLOCK_SIZE.
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
                                perm_stats_t *pperm_stats,
                                float *pES_perm_vec,
                                const float *pES,
                                const ES_bounds_t *pES_bounds,
                                const int n_sizes,
                                const int block_size) {
  for (int s = 0; s < n_sizes; ++s) {
    const int block_start_perm = s * block_size;

    // Partition permutation ES by sign and update the sums ====
    int first_positive_perm = block_start_perm;
    float sum_negative_perm = 0.0f;
    float sum_positive_perm = 0.0f;

    for (int k = 0; k < block_size; ++k) {
      const float swap = pES_perm_vec[block_start_perm + k];

      if (swap < 0.0f) {
        pES_perm_vec[block_start_perm + k] = pES_perm_vec[first_positive_perm];
        pES_perm_vec[first_positive_perm] = swap;

        sum_negative_perm -= swap;
        ++first_positive_perm;
      } else {
        sum_positive_perm += swap;
      }
    }

    pperm_stats[s].n_negative += (first_positive_perm - block_start_perm);
    pperm_stats[s].sum_negative += sum_negative_perm;
    pperm_stats[s].sum_positive += sum_positive_perm;

    // Comparisons ====
    int i_neg = pES_bounds[s].start;
    int i_pos = pES_bounds[s].first_positive;
    int remainder_neg = i_pos - i_neg;
    int remainder_pos = pES_bounds[s].end - i_pos;
    const int perm_end = block_start_perm + block_size;

    // Negative ----
    while (remainder_neg >= 8) {
      for (int b = block_start_perm; b < first_positive_perm; ++b) {
        for (int k = 0; k < 8; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[b] <= pES[i_neg + k];
        }
      }
      i_neg += 8;
      remainder_neg -= 8;
    }

    if (remainder_neg >= 4) {
      for (int b = block_start_perm; b < first_positive_perm; ++b) {
        for (int k = 0; k < 4; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[b] <= pES[i_neg + k];
        }
      }
      i_neg += 4;
      remainder_neg -= 4;
    }

    if (remainder_neg >= 2) {
      for (int b = block_start_perm; b < first_positive_perm; ++b) {
        for (int k = 0; k < 2; ++k) {
          pn_as_extreme[i_neg + k] += pES_perm_vec[b] <= pES[i_neg + k];
        }
      }
      i_neg += 2;
      remainder_neg -= 2;
    }

    if (remainder_neg == 1) {
      for (int b = block_start_perm; b < first_positive_perm; ++b) {
        pn_as_extreme[i_neg] += pES_perm_vec[b] <= pES[i_neg];
      }
    }

    // Positive ----
    while (remainder_pos >= 8) {
      for (int b = first_positive_perm; b < perm_end; ++b) {
        for (int k = 0; k < 8; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[b] >= pES[i_pos + k];
        }
      }
      i_pos += 8;
      remainder_pos -= 8;
    }

    if (remainder_pos >= 4) {
      for (int b = first_positive_perm; b < perm_end; ++b) {
        for (int k = 0; k < 4; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[b] >= pES[i_pos + k];
        }
      }
      i_pos += 4;
      remainder_pos -= 4;
    }

    if (remainder_pos >= 2) {
      for (int b = first_positive_perm; b < perm_end; ++b) {
        for (int k = 0; k < 2; ++k) {
          pn_as_extreme[i_pos + k] += pES_perm_vec[b] >= pES[i_pos + k];
        }
      }
      i_pos += 2;
      remainder_pos -= 2;
    }

    if (remainder_pos == 1) {
      for (int b = first_positive_perm; b < perm_end; ++b) {
        pn_as_extreme[i_pos] += pES_perm_vec[b] >= pES[i_pos];
      }
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

  std::vector<gene_data_t> gene_data(n_genes);
  gene_data_t *pgene_data = &gene_data[0];
  for (int i = 0; i < n_genes; ++i) {
    pgene_data[i].r = (float) pr_dbl[i];
    pgene_data[i].y = (float) py_dbl[i];
  }

  std::vector<perm_stats_t> perm_stats(n_sizes);
  perm_stats_t *pperm_stats = &perm_stats[0];

  std::vector<ES_bounds_t> ES_bounds(n_sizes);
  ES_bounds_t *pES_bounds = &ES_bounds[0];
  pES_bounds[0].start = 0;
  pES_bounds[0].end = pES_end[0];
  for (int i = 1; i < n_sizes; ++i) {
    pES_bounds[i].start = pES_end[i - 1];
    pES_bounds[i].end = pES_end[i];
  }
  update_first_positive(pES_bounds, n_sizes, pES);

  std::vector<float> inv_w(n_sizes);
  float *pinv_w = &inv_w[0];
  for (int i = 0; i < n_sizes; ++i) {
    pinv_w[i] = (float) (1.0 / punique_w[i]);
  }

  std::vector<float> ES_perm_vec(BLOCK_SIZE * n_sizes);
  float *pES_perm_vec = &ES_perm_vec[0];

  dqrng::dqset_seed(seed);
  const int partial_block_size = nperm % BLOCK_SIZE;

  if (partial_block_size != 0) {
    calc_ES_perm_internal(
      pES_perm_vec,
      n_genes,
      n_sizes,
      partial_block_size,
      pgene_data,
      max_size,
      sum_ranks,
      punique_m,
      pinv_w
    );

    update_n_as_extreme(
      pn_as_extreme,
      pperm_stats,
      pES_perm_vec,
      pES,
      pES_bounds,
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
      pgene_data,
      max_size,
      sum_ranks,
      punique_m,
      pinv_w
    );

    update_n_as_extreme(
      pn_as_extreme,
      pperm_stats,
      pES_perm_vec,
      pES,
      pES_bounds,
      n_sizes,
      BLOCK_SIZE
    );
  }

  update_output_vectors(
    pn_same_sign,
    psum_ES_perm,
    n_sizes,
    nperm,
    pES_bounds,
    pperm_stats
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
                                      const gene_data_t *pgene_data,
                                      const int max_size,
                                      const float sum_ranks,
                                      const int *punique_m_down,
                                      const int *punique_m_up,
                                      const float *pinv_w_up,
                                      const float *pinv_w_down,
                                      const pair_map_t *ppair_map) {
  for (int b = 0; b < block_size; ++b) {
    const Rcpp::IntegerVector random_indices = dqrng::dqsample_int(
      n_genes,
      max_size
    );
    const int *prandom_indices = random_indices.begin();

    // Up-regulated genes ====
    int start = 0;
    float sum_r = 0.0f;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int i_up = 0; i_up < n_sizes_up; ++i_up) {
      const int end = punique_m_up[i_up];
      const bool empty_set = (start == end);

      while (start < end) {
        const int rand_idx = prandom_indices[start++];
        sum_r += pgene_data[rand_idx].r;
        sum_y += pgene_data[rand_idx].y;
        sum_ry += pgene_data[rand_idx].r * pgene_data[rand_idx].y;
      }

      pES_perm_up[i_up] = empty_set ? 0.0f :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_w_up[i_up];
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
      // end of random_indices to avoid selecting the same values that were used
      // for the up-regulated permutation ES.
      while (start > end) {
        const int rand_idx = prandom_indices[start--];
        sum_r += pgene_data[rand_idx].r;
        sum_y += pgene_data[rand_idx].y;
        sum_ry += pgene_data[rand_idx].r * pgene_data[rand_idx].y;
      }

      pES_perm_down[i_down] = empty_set ? 0.0f :
        (sum_ry / sum_y) + (sum_r - sum_ranks) * pinv_w_down[i_down];
    }

    // ES = ES_up - ES_down for each unique pair
    int perm_idx = b;
    for (int s = 0; s < n_pairs; ++s, perm_idx += block_size) {
      pES_perm_vec[perm_idx] = pES_perm_up[ppair_map[s].up] - pES_perm_down[ppair_map[s].down];
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

  std::vector<gene_data_t> gene_data(n_genes);
  gene_data_t *pgene_data = &gene_data[0];
  for (int i = 0; i < n_genes; ++i) {
    pgene_data[i].r = (float) pr_dbl[i];
    pgene_data[i].y = (float) py_dbl[i];
  }

  std::vector<pair_map_t> pair_map(n_pairs);
  pair_map_t *ppair_map = &pair_map[0];
  // Convert from 1-based to 0-based indices
  for (int i = 0; i < n_pairs; ++i) {
    ppair_map[i].up = pmap_unique_to_pairs_up[i] - 1;
    ppair_map[i].down = pmap_unique_to_pairs_down[i] - 1;
  }

  std::vector<perm_stats_t> perm_stats(n_pairs);
  perm_stats_t *pperm_stats = &perm_stats[0];

  std::vector<ES_bounds_t> ES_bounds(n_pairs);
  ES_bounds_t *pES_bounds = &ES_bounds[0];
  pES_bounds[0].start = 0;
  pES_bounds[0].end = pES_end[0];
  for (int i = 1; i < n_pairs; ++i) {
    pES_bounds[i].start = pES_end[i - 1];
    pES_bounds[i].end = pES_end[i];
  }
  update_first_positive(pES_bounds, n_pairs, pES);

  // Update end positions of down-regulated genes to avoid overlap with
  // up-regulated genes
  for (int i = 0; i < n_sizes_down; ++i) {
    punique_m_down[i] = (max_size - 1) - punique_m_down[i];
  }

  std::vector<float> inv_w_up(n_sizes_up);
  float *pinv_w_up = &inv_w_up[0];
  for (int i = 0; i < n_sizes_up; ++i) {
    pinv_w_up[i] = (float) (1.0 / punique_w_up[i]);
  }

  std::vector<float> inv_w_down(n_sizes_down);
  float *pinv_w_down = &inv_w_down[0];
  for (int i = 0; i < n_sizes_down; ++i) {
    pinv_w_down[i] = (float) (1.0 / punique_w_down[i]);
  }

  std::vector<float> ES_perm_up(n_sizes_up);
  std::vector<float> ES_perm_down(n_sizes_down);
  std::vector<float> ES_perm_vec(BLOCK_SIZE * n_pairs);
  float *pES_perm_up = &ES_perm_up[0];
  float *pES_perm_down = &ES_perm_down[0];
  float *pES_perm_vec = &ES_perm_vec[0];

  dqrng::dqset_seed(seed);
  const int partial_block_size = nperm % BLOCK_SIZE;

  if (partial_block_size != 0) {
    calc_ES_perm_dir_internal(
      pES_perm_vec,
      pES_perm_down,
      pES_perm_up,
      n_genes,
      n_sizes_down,
      n_sizes_up,
      n_pairs,
      partial_block_size,
      pgene_data,
      max_size,
      sum_ranks,
      punique_m_down,
      punique_m_up,
      pinv_w_up,
      pinv_w_down,
      ppair_map
    );

    update_n_as_extreme(
      pn_as_extreme,
      pperm_stats,
      pES_perm_vec,
      pES,
      pES_bounds,
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
      pgene_data,
      max_size,
      sum_ranks,
      punique_m_down,
      punique_m_up,
      pinv_w_up,
      pinv_w_down,
      ppair_map
    );

    update_n_as_extreme(
      pn_as_extreme,
      pperm_stats,
      pES_perm_vec,
      pES,
      pES_bounds,
      n_pairs,
      BLOCK_SIZE
    );
  }

  update_output_vectors(
    pn_same_sign,
    psum_ES_perm,
    n_pairs,
    nperm,
    pES_bounds,
    pperm_stats
  );
}
