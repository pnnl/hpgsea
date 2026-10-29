#include <Rcpp.h>
// [[Rcpp::depends(dqrng, BH)]]
#include <dqrng_generator.h>
#include <dqrng_sample.h>
#include <dqrng.h>

#ifdef _OPENMP
  #include <omp.h>
  #define GET_THREAD_ID() omp_get_thread_num()
  #define OMP_MAX_THREADS omp_get_max_threads()
#else
  #define GET_THREAD_ID() 0
  #define OMP_MAX_THREADS 1
#endif

// [[Rcpp::plugins(openmp)]]

#define BLOCK_SIZE 200 // permutations per block
#define TOWER_SIZE 10000 // permutations per tower
#define TOWERS_PER_HYPERTOWER 5

#define USE_SIMD

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

dqrng::rng64_t rng = dqrng::generator();

// Gene-level values
typedef struct {
  float r;
  float y;
} GeneData;

// For a given set size, defines the index of the first ES, the index of the
// first positive ES, and the index of the last ES + 1.
typedef struct {
  int start;
  int first_positive;
  int end;
} ESBounds;

// For a given set size, count the number of negative permutation ES and compute
// the absolute sums of the negative and positive ES, separately.
typedef struct {
  int n_negative;
  float sum_negative;
  float sum_positive;
} PermStats;

typedef struct {
  int up;
  int down;
} PairMap;


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
void update_first_positive(ESBounds *pES_bounds,
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
                                  const ESBounds *pES_bounds,
                                  const PermStats *pperm_stats)
{
  for (int s = 0; s < n_sizes; ++s) {
    const int start = pES_bounds[s].start;
    const int first_positive = pES_bounds[s].first_positive;
    const int end = pES_bounds[s].end;

    const int n_negative = pperm_stats[s].n_negative;
    const int n_positive = nperm - n_negative;
    const double sum_negative = (double) -pperm_stats[s].sum_negative;
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


inline void calc_ES_perm_internal(dqrng::random_64bit_generator &_rng,
                                  float *pES_perm_vec,
                                  const int n_genes,
                                  const int n_sizes,
                                  const int block_size,
                                  const GeneData *pgene_data,
                                  const int max_size,
                                  const float sum_ranks,
                                  const int *punique_m,
                                  const float *pinv_w) {
  for (int b = 0; b < block_size; ++b) {
    std::vector<int> random_indices =
      dqrng::sample::sample<std::vector<int>, uint32_t>(
        _rng, uint32_t(n_genes), uint32_t(max_size), false, 0
      );
    const int *prandom_indices = &random_indices[0];

    int perm_idx = b;
    int start = 0;
    float delta_ranks = -sum_ranks;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int s = 0; s < n_sizes; ++s, perm_idx += block_size) {
      const int end = punique_m[s];

      while (start < end) {
        const GeneData gene = pgene_data[prandom_indices[start++]];
        delta_ranks += gene.r;
        sum_y += gene.y;
        sum_ry += gene.r * gene.y;
      }

      pES_perm_vec[perm_idx] = (sum_ry / sum_y) + delta_ranks * pinv_w[s];
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

#ifdef USE_SIMD

inline void update_n_as_extreme(int *pn_as_extreme,
                                PermStats *pperm_stats,
                                float *pES_perm_vec,
                                const float *pES,
                                const ESBounds *pES_bounds,
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

        sum_negative_perm += swap;
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

#else

inline void update_n_as_extreme(int *pn_as_extreme,
                                PermStats *pperm_stats,
                                const float *pES_perm_vec,
                                const float *pES,
                                const ESBounds *pES_bounds,
                                const int n_sizes,
                                const int block_size) {
  for (int s = 0; s < n_sizes; ++s) {
    const int block_start_perm = s * block_size;

    // Partition permutation ES by sign and update the sums ====
    int n_negative = block_start_perm;
    float sum_negative_perm = 0.0f;
    float sum_positive_perm = 0.0f;

    for (int k = 0; k < block_size; ++k) {
      const float swap = pES_perm_vec[block_start_perm + k];

      if (swap < 0.0f) {
        sum_negative_perm += swap;
        ++n_negative;
      } else {
        sum_positive_perm += swap;
      }
    }

    pperm_stats[s].n_negative += (n_negative - block_start_perm);
    pperm_stats[s].sum_negative += sum_negative_perm;
    pperm_stats[s].sum_positive += sum_positive_perm;

    // Comparisons ====
    const int start_s = pES_bounds[s].start;
    const int end_s = pES_bounds[s].end;

    for (int i = start_s; i < end_s; ++i) {
      int total = 0;
      const float ES_i = pES[i];

      if (ES_i < 0.0f) {
        for (int j = 0; j < block_size; ++j) {
          total += pES_perm_vec[block_start_perm + j] <= ES_i;
        }
      } else {
        for (int j = 0; j < block_size; ++j) {
          total += pES_perm_vec[block_start_perm + j] >= ES_i;
        }
      }

      pn_as_extreme[i] += total;
    }

  }
}

#endif


uint64_t get_tower_stream(int tower_id) {
  return (uint64_t) (tower_id + 1);
}


void combine_thread_results(const int nthreads,
                            const int n_towers,
                            const int n_sets,
                            const int n_sizes,
                            PermStats *pperm_stats,
                            int *pn_as_extreme_threads,
                            int *pn_as_extreme) {
  for (int j = 1; j < nthreads; ++j) {
    const int offset1 = n_sizes * j;

    for (int i = 0; i < n_sizes; ++i) {
      pperm_stats[i].n_negative += pperm_stats[i + offset1].n_negative;
      pperm_stats[i].sum_negative += pperm_stats[i + offset1].sum_negative;
      pperm_stats[i].sum_positive += pperm_stats[i + offset1].sum_positive;
    }
  }

  for (int j = 1; j < nthreads; ++j) {
    const int offset2 = n_sets * j;

    for (int i = 0; i < n_sets; ++i) {
      pn_as_extreme_threads[i] += pn_as_extreme_threads[i + offset2];
    }
  }

  if (n_towers != 0) {
    for (int i = 0; i < n_sets; ++i) {
      pn_as_extreme[i] += pn_as_extreme_threads[i];
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
                  const int nperm,
                  const SEXP ES_dbl,
                  const SEXP ES_end,
                  const SEXP y_dbl,
                  const SEXP r_dbl,
                  const int max_size,
                  const SEXP Rsum_ranks,
                  const SEXP unique_m,
                  const SEXP unique_w,
                  const SEXP Rnthreads)
{
  const int n_towers = nperm / TOWER_SIZE;
  int nthreads = MIN(OMP_MAX_THREADS, INTEGER(Rnthreads)[0]);
  nthreads = n_towers < nthreads ? MAX(n_towers, 1) : nthreads;
  nthreads = (nthreads == 0) ? MAX(1, MIN(n_towers, OMP_MAX_THREADS)) : nthreads;

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

  std::vector<GeneData> gene_data(n_genes);
  GeneData *pgene_data = &gene_data[0];
  for (int i = 0; i < n_genes; ++i) {
    pgene_data[i].r = (float) pr_dbl[i];
    pgene_data[i].y = (float) py_dbl[i];
  }

  std::vector<PermStats> perm_stats(n_sizes * nthreads);
  PermStats *pperm_stats = &perm_stats[0];

  std::vector<ESBounds> ES_bounds(n_sizes);
  ESBounds *pES_bounds = &ES_bounds[0];
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

  const int partial_tower_size = (nperm < TOWER_SIZE) ? nperm :
    (nperm % TOWER_SIZE);

  std::vector<float> ES_perm_vec(n_sizes * BLOCK_SIZE * nthreads);
  float *pES_perm_vec = &ES_perm_vec[0];

  std::vector<int> n_as_extreme_threads(n_sets * nthreads);
  int *pn_as_extreme_threads = &n_as_extreme_threads[0];

  uint64_t _seed = dqrng::get_seed_from_r();
  rng->seed(_seed, (uint64_t) 0);

  // The first tower may be partially filled. The individual blocks in the tower
  // may be complete, or the first block may be a partial one.
  if (partial_tower_size != 0) {
    const int partial_block_size = (partial_tower_size < BLOCK_SIZE) ?
      partial_tower_size : (partial_tower_size % BLOCK_SIZE);

    if (partial_block_size != 0) {
      calc_ES_perm_internal(
        *rng,
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

    Rcpp::checkUserInterrupt();

    const int n_blocks = partial_tower_size / BLOCK_SIZE;

    for (int i = 0; i < n_blocks; ++i) {
      calc_ES_perm_internal(
        *rng,
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
  }

  const int n_blocks_per_tower = TOWER_SIZE / BLOCK_SIZE;

  int tower_id = 0;
  int towers_remaining = n_towers;

  // I create hypertowers (collections of towers) because
  // Rcpp::checkUserInterrupt() must be called outside of the OpenMP for loop.
  const int hypertower_size = nthreads * TOWERS_PER_HYPERTOWER;

  // Complete hypertowers
  while (towers_remaining >= hypertower_size) {
    Rcpp::checkUserInterrupt();

    #pragma omp parallel for num_threads(nthreads)
    for (int _tower_id = 0; _tower_id < hypertower_size; ++_tower_id) {
      const int thread_id = GET_THREAD_ID();

      uint64_t tower_stream = get_tower_stream(tower_id + _tower_id);
      auto rng_tower = rng->clone(tower_stream);

      // Pointers to thread-safe chunks of data
      float *pES_perm_vec_tower = pES_perm_vec + thread_id * BLOCK_SIZE * n_sizes;
      PermStats *pperm_stats_tower = pperm_stats + thread_id * n_sizes;
      int *pn_as_extreme_tower = pn_as_extreme_threads + thread_id * n_sets;

      for (int i = 0; i < n_blocks_per_tower; ++i) {
        calc_ES_perm_internal(
          *rng_tower,
          pES_perm_vec_tower,
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
          pn_as_extreme_tower,
          pperm_stats_tower,
          pES_perm_vec_tower,
          pES,
          pES_bounds,
          n_sizes,
          BLOCK_SIZE
        );
      }
    }

    tower_id += hypertower_size;
    towers_remaining -= hypertower_size;
  }

  Rcpp::checkUserInterrupt();

  // Partial hypertower
  #pragma omp parallel for num_threads(nthreads)
  for (int _tower_id = 0; _tower_id < towers_remaining; ++_tower_id) {
    const int thread_id = GET_THREAD_ID();

    uint64_t tower_stream = get_tower_stream(tower_id + _tower_id);
    auto rng_tower = rng->clone(tower_stream);

    // Pointers to thread-safe chunks of data
    float *pES_perm_vec_tower = pES_perm_vec + thread_id * BLOCK_SIZE * n_sizes;
    PermStats *pperm_stats_tower = pperm_stats + thread_id * n_sizes;
    int *pn_as_extreme_tower = pn_as_extreme_threads + thread_id * n_sets;

    for (int i = 0; i < n_blocks_per_tower; ++i) {
      calc_ES_perm_internal(
        *rng_tower,
        pES_perm_vec_tower,
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
        pn_as_extreme_tower,
        pperm_stats_tower,
        pES_perm_vec_tower,
        pES,
        pES_bounds,
        n_sizes,
        BLOCK_SIZE
      );
    }
  }

  combine_thread_results(
    nthreads,
    n_towers,
    n_sets,
    n_sizes,
    pperm_stats,
    pn_as_extreme_threads,
    pn_as_extreme
  );

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
                                      const GeneData *pgene_data,
                                      const int max_size,
                                      const float sum_ranks,
                                      const int *punique_m_down,
                                      const int *punique_m_up,
                                      const float *pinv_w_up,
                                      const float *pinv_w_down,
                                      const PairMap *ppair_map) {
  int perm_idx;

  for (int b = 0; b < block_size; ++b) {
    std::vector<int> random_indices =
      dqrng::sample::sample<std::vector<int>, uint32_t>(
          *rng, uint32_t(n_genes), uint32_t(max_size), false, 0
      );
    const int *prandom_indices = &random_indices[0];

    // Up-regulated genes ====
    perm_idx = b;
    int start = 0;
    float delta_ranks = -sum_ranks;
    float sum_y = 0.0f;
    float sum_ry = 0.0f;

    for (int i_up = 0; i_up < n_sizes_up; ++i_up, perm_idx += block_size) {
      const int end = punique_m_up[i_up];
      const bool empty_set = (start == end);

      while (start < end) {
        const GeneData gene = pgene_data[prandom_indices[start++]];
        delta_ranks += gene.r;
        sum_y += gene.y;
        sum_ry += gene.r * gene.y;
      }

      pES_perm_up[perm_idx] = empty_set ? 0.0f :
        (sum_ry / sum_y) + delta_ranks * pinv_w_up[i_up];
    }

    // Down-regulated genes ====
    perm_idx = b;
    start = max_size - 1;
    delta_ranks = -sum_ranks;
    sum_y = 0.0f;
    sum_ry = 0.0f;

    for (int i_down = 0; i_down < n_sizes_down; ++i_down, perm_idx += block_size) {
      const int end = punique_m_down[i_down];
      const bool empty_set = (start == end);

      // Indices for the down-regulated portions of the gene sets start at the
      // end of random_indices to avoid selecting the same values that were used
      // for the up-regulated permutation ES.
      while (start > end) {
        const GeneData gene = pgene_data[prandom_indices[start--]];
        delta_ranks += gene.r;
        sum_y += gene.y;
        sum_ry += gene.r * gene.y;
      }

      pES_perm_down[perm_idx] = empty_set ? 0.0f :
        (sum_ry / sum_y) + delta_ranks * pinv_w_down[i_down];
    }
  }

  perm_idx = 0;

  // ES = ES_up - ES_down for each unique pair
  for (int s = 0; s < n_pairs; ++s) {
    const PairMap pair = ppair_map[s];

    for (int b = 0; b < block_size; ++b) {
      pES_perm_vec[perm_idx++] =
        pES_perm_up[pair.up * block_size + b] -
        pES_perm_down[pair.down * block_size + b];
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

  std::vector<GeneData> gene_data(n_genes);
  GeneData *pgene_data = &gene_data[0];
  for (int i = 0; i < n_genes; ++i) {
    pgene_data[i].r = (float) pr_dbl[i];
    pgene_data[i].y = (float) py_dbl[i];
  }

  std::vector<PairMap> pair_map(n_pairs);
  PairMap *ppair_map = &pair_map[0];
  // Convert from 1-based to 0-based indices
  for (int i = 0; i < n_pairs; ++i) {
    ppair_map[i].up = pmap_unique_to_pairs_up[i] - 1;
    ppair_map[i].down = pmap_unique_to_pairs_down[i] - 1;
  }

  std::vector<PermStats> perm_stats(n_pairs);
  PermStats *pperm_stats = &perm_stats[0];

  std::vector<ESBounds> ES_bounds(n_pairs);
  ESBounds *pES_bounds = &ES_bounds[0];
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

  std::vector<float> ES_perm_up(BLOCK_SIZE * n_sizes_up);
  std::vector<float> ES_perm_down(BLOCK_SIZE * n_sizes_down);
  std::vector<float> ES_perm_vec(BLOCK_SIZE * n_pairs);
  float *pES_perm_up = &ES_perm_up[0];
  float *pES_perm_down = &ES_perm_down[0];
  float *pES_perm_vec = &ES_perm_vec[0];

  const int partial_block_size = nperm % BLOCK_SIZE;

  uint64_t _seed = dqrng::get_seed_from_r();
  rng->seed(_seed);

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
