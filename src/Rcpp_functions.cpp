// [[Rcpp::depends(RcppArmadillo, dqrng)]]
#include <RcppArmadillo.h>
#include <dqrng.h>

using namespace Rcpp;


// [[Rcpp::export(.Cpp_unsafe_sparseMatrix)]]
SEXP unsafe_sparseMatrix(IntegerVector& i,
                         const IntegerVector& j,
                         const IntegerVector dims,
                         const List& dimnames) {
  // j is sorted and i is sorted within j. Duplicate (i, j) pairs not allowed.
  const int N = i.size();
  const int NCOL = dims[1];

  // Row indices need to be 0-based
  for (int k = 0; k < N; ++k) {
    --i[k];
  }

  // Column pointers
  IntegerVector p(NCOL + 1, 0);

  for (int k = 0; k < N; ++k) {
    const int col = j[k];
    ++p[col]; // column sums
  }

  for (int col = 0; col < NCOL; ++col) {
    p[col + 1] += p[col]; // cumulative sum
  }

  // Create a dgCMatrix
  S4 out("dgCMatrix");
  out.slot("i") = i;
  out.slot("p") = p;
  out.slot("x") = NumericVector(N, 1.0);
  out.slot("Dim") = dims;
  out.slot("Dimnames") = dimnames;

  return out;
}


//' @title Dense-Sparse Matrix Multiplication
//'
//' @description Multiply a dense matrix by a sparse matrix.
//'
//' @param X dense matrix.
//' @param Y sparse matrix.
//'
//' @returns The product of \code{X} and \code{Y}: a dense matrix.
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//' @noRd

// [[Rcpp::export(.Cpp_matmult_sparse)]]
arma::mat matmult_sparse(const arma::mat& X,
                         const arma::sp_mat& Y)
{
  return X * Y;
}


//' @title Core of the .calcES R function
//'
//' @param Y absolute values of the matrix \code{t(X)} raised to the power of
//'   \code{alpha}. Missing values are then imputed with 0.
//' @param R matrix of ranks of the values in each row of \code{X}. Missing
//'   values in \code{X} are assigned a rank of \code{NA}, which are then
//'   imputed with 0.
//' @param sumRanks integer vector; the sums of the ranks in each sample. Equal
//'   to \code{rowSums(R)}.
//' @param A sparse incidence matrix with single genes as rows and gene sets as
//'   columns. A value of 1 indicates that the gene is an element of the set,
//'   while a value of 0 indicates otherwise.
//' @param M matrix with samples as rows and gene sets as columns, where each
//'   entry is the number of genes with nonmissing values in each set.
//' @param W matrix with the same dimensions as \code{M} where each entry is the
//'   number of genes with nonmissing values \emph{not} in each set.
//'
//' @returns A matrix of real-valued enrichment scores with samples as rows and
//'   gene sets as columns. May contain missing values if the corresponding
//'   entry of \code{M} is less than 2.
//'
//' @author Tyler Sagendorf
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//'   Sanderson, C., & Curtin, R. (2019). Practical Sparse Matrices in C++ with
//'   Hybrid Storage and Template-Based Expression Optimisation. Mathematical
//'   and Computational Applications, 24(3), 70. \url{
//'   https://doi.org/10.3390/mca24030070}
//'
//' @noRd

// [[Rcpp::export(.Cpp_calcES)]]
arma::mat calcES(const int min_size,
                 const arma::mat& Y,
                 const arma::mat& R,
                 const arma::colvec& sumRanks,
                 const arma::sp_mat& A,
                 const arma::mat& M,
                 const arma::mat& W)
{
  // % is Hadamard product, * is dot product, / is Hadamard division
  arma::mat RA = R * A;
  RA.each_col() -= sumRanks;

  arma::mat ES = ((R % Y) * A) / (Y * A) + (RA / W);

  // If the set has fewer than min_size elements with nonmissing values, the ES
  // will be 0. Only applies to directional gene sets.
  arma::uvec indices = arma::find(M < min_size);
  ES(indices).zeros();

  return ES;
}


// Permutation Tests ===========================================================

//' @title Get the Index of the First Positive ES for Every Unique Gene Set Size
//'
//' @description Perform binary searches within sections of a numeric vector to
//'   determine the index of the first positive ES.
//'
//' @param n_sizes integer; the number of unique gene set sizes.
//' @param ES numeric vector of enrichment scores, sorted in ascending order by
//'   gene set size and then by the values of the ES.
//' @param ES_start integer vector with length equal to \code{n_sizes}. Stores
//'   the index of the first ES for every unique gene set size.
//' @param ES_end integer vector with length equal to \code{n_sizes}. Each
//'   element is 1 more than the index of the last ES for every unique gene set
//'   size.
//'
//' @returns An integer vector containing the index of the first positive ES for
//'   every unique gene set size. Each element is greater than or equal to the
//'   corresponding element of \code{ES_start}. If there are no positive ES for
//'   a particular set size, the index will be the corresponding element of
//'   \code{ES_end}.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
IntegerVector get_ES_pos_idx(const int n_sizes,
                             NumericVector ES,
                             IntegerVector ES_start,
                             IntegerVector ES_end)
{
  // Indices of first positive ES by unique set size
  IntegerVector ES_pos_idx(n_sizes, 0);

  // zero_idx is the minimum index for the ES from gene sets of a particular
  // size. Prevents checking ES from smaller gene sets.
  int zero_idx, low, mid, high;

  for (int i = 0; i < n_sizes; ++i) {
    low = ES_start[i];
    zero_idx = low;
    high = ES_end[i] - 1;
    ES_pos_idx[i] = ES_end[i];

    while (low <= high) {
      mid = (low + high) >> 1;

      if (ES[mid] >= 0.0) {
        if (mid == zero_idx || ES[mid - 1] < 0.0) {
          ES_pos_idx[i] = mid;
          break; // break out of while loop
        } else {
          high = mid - 1;
        }
      } else {
        low = mid + 1;
      }
    }
  }

  return ES_pos_idx;
}


// From src/C_functions.c
extern "C" {
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
                      SEXP inv_L2_w);
}


//' @title Update n_same_sign and sum_ES_perm vectors
//'
//' @description Map the result vectors with length equal to the number of
//'   unique gene set sizes to the result vectors with length equal to the total
//'   number of gene sets.
//'
//' @param n_same_sign integer vector of zeros with length equal to the number
//'   of enrichment scores. Stores the number of permutation ES with the same
//'   sign as each true ES.
//' @param sum_ES_perm numeric vector of zeros with length equal to the number
//'   of enrichment scores. Stores the absolute sum of the permutation ES with
//'   the same sign as each true ES.
//' @param n_sizes integer; number of unique gene set sizes.
//' @param nperm integer; total number of permutations.
//' @param ES_start integer vector with length equal to \code{n_sizes}. Stores
//'   the index of the first ES for every unique gene set size.
//' @param ES_end integer vector with length equal to \code{n_sizes}. Each
//'   element is 1 more than the index of the last ES for every unique gene set
//'   size.
//' @param ES_pos_idx integer vector; the output of \code{get_ES_pos_idx()}.
//' @param n_perm_neg integer vector; contains the number of negative
//'   permutation ES for every unique gene set size.
//' @param sum_perm_pos numeric vector; the sum of the positive permutation ES
//'   for every unique gene set size.
//' @param sum_perm_neg numeric vector; the absolute sum of the negative
//'   permutation ES for every unique gene set size.
//'
//' @returns Nothing. The vectors \code{n_same_sign} and \code{sum_ES_perm} are
//'   modified in place.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
void update_results(IntegerVector n_same_sign,
                    NumericVector sum_ES_perm,
                    const int n_sizes,
                    const int nperm,
                    IntegerVector ES_start,
                    IntegerVector ES_end,
                    IntegerVector ES_pos_idx,
                    IntegerVector n_perm_neg,
                    NumericVector sum_perm_pos,
                    NumericVector sum_perm_neg)
{
  int ES_start_i, ES_end_i, ES_pos_idx_i, n_perm_neg_i, n_perm_pos_i;
  double sum_perm_neg_i, sum_perm_pos_i;

  for (int i = 0; i < n_sizes; ++i) {
    ES_start_i = ES_start[i];
    ES_end_i = ES_end[i];
    ES_pos_idx_i = ES_pos_idx[i];

    n_perm_neg_i = n_perm_neg[i];
    n_perm_pos_i = nperm - n_perm_neg_i;

    sum_perm_neg_i = sum_perm_neg[i];
    sum_perm_pos_i = sum_perm_pos[i];

    for (int j = ES_start_i; j < ES_pos_idx_i; ++j) {
      n_same_sign[j] = n_perm_neg_i;
      sum_ES_perm[j] = sum_perm_neg_i;
    }

    for (int j = ES_pos_idx_i; j < ES_end_i; ++j) {
      n_same_sign[j] = n_perm_pos_i;
      sum_ES_perm[j] = sum_perm_pos_i;
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
//' @param batch_size integer; the number of permutations run as a single batch.
//' @param ES numeric vector of enrichment scores, sorted in ascending order by
//'   gene set size and then by the values of the ES.
//' @param ES_end integer vector with length equal to \code{n_sizes}. Each
//'   element is 1 more than the index of the last ES for every unique gene set
//'   size.
//' @param y the absolute values of the gene level statistics raised to some
//'   non-negative power \code{alpha}.
//' @param r the ranks of the gene-level statistics.
//' @param max_set_size integer; the size of the largest gene set.
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
void calc_ES_perm(IntegerVector& n_same_sign,
                  IntegerVector& n_as_extreme,
                  NumericVector& sum_ES_perm,
                  const Nullable<IntegerVector> seed,
                  const int nperm,
                  int batch_size,
                  const NumericVector& ES,
                  const IntegerVector& ES_end,
                  const NumericVector& y,
                  const NumericVector& r,
                  const int max_set_size,
                  const double sum_ranks,
                  const IntegerVector& L2_m,
                  const IntegerVector& L2_w)
{
  const int n_genes = y.length(); // number of gene-level statistics
  const int n_sizes = L2_m.length(); // number of unique gene set sizes

  // Vectors updated by C_calc_ES_perm (n_as_extreme is also updated)
  IntegerVector n_perm_neg(n_sizes, 0);
  NumericVector sum_perm_pos(n_sizes, 0.0);
  NumericVector sum_perm_neg(n_sizes, 0.0);

  // The NumericVector version of L2_w is needed later to divide quantities when
  // calculating permutation ES anyway, so may as well take its inverse here so
  // that multiplication can be used instead of division, since multiplication
  // takes fewer CPU cycles.
  NumericVector inv_L2_w(n_sizes);

  for (int i = 0; i < n_sizes; ++i) {
    inv_L2_w[i] = 1.0 / (double)(L2_w[i]);
  }

  IntegerVector ES_start(n_sizes, 0);

  for (int i = 1; i < n_sizes; ++i) {
    ES_start[i] = ES_end[i - 1];
  }

  // Index of first positive ES by unique set size
  IntegerVector ES_pos_idx = get_ES_pos_idx(n_sizes, ES, ES_start, ES_end);

  // Starting positions for segmented sums
  IntegerVector start_vec(n_sizes, 0);

  for (int i = 1; i < n_sizes; ++i) {
    start_vec[i] = L2_m[i - 1];
  }

  int n_batches = std::ceil((double)nperm / (double)batch_size);

  // The first batch may be smaller than the rest if nperm is not an integer
  // multiple of n_batches
  const int leftover_perm = nperm % n_batches;

  dqrng::dqset_seed(seed);
  IntegerMatrix perm_idx(max_set_size, batch_size);

  if (leftover_perm > 0) {
    // Each column is a random sample of size max_set_size from integers 0 to
    // n_genes - 1. Used to select pairs of elements from the vectors y and r.
    for (int j = 0; j < leftover_perm; ++j) {
      perm_idx.column(j) = dqrng::dqsample_int(n_genes, max_set_size);
    }

    C_calc_ES_perm(
      n_as_extreme,
      n_perm_neg,
      sum_perm_pos,
      sum_perm_neg,
      ES,
      ES_start,
      ES_end,
      ES_pos_idx,
      y,
      r,
      perm_idx,
      max_set_size,
      leftover_perm, // batch_size
      sum_ranks,
      start_vec,
      L2_m, // end_vec
      inv_L2_w
    );

    --n_batches;
  }

  // Loop over remaining permutations (integer multiple of n_batches)
  for (int batch = 0; batch < n_batches; ++batch) {
    for (int j = 0; j < batch_size; ++j) {
      perm_idx.column(j) = dqrng::dqsample_int(n_genes, max_set_size);
    }

    C_calc_ES_perm(
      n_as_extreme,
      n_perm_neg,
      sum_perm_pos,
      sum_perm_neg,
      ES,
      ES_start,
      ES_end,
      ES_pos_idx,
      y,
      r,
      perm_idx,
      max_set_size,
      batch_size,
      sum_ranks,
      start_vec,
      L2_m,
      inv_L2_w
    );
  }

  // Update n_same_sign and sum_ES_perm. n_as_extreme was updated by
  // C_calc_ES_perm.
  update_results(n_same_sign,
                 sum_ES_perm,
                 n_sizes,
                 nperm,
                 ES_start,
                 ES_end,
                 ES_pos_idx,
                 n_perm_neg,
                 sum_perm_pos,
                 sum_perm_neg);
}


// Directional Gene Set Permutation Tests ======================================

// From src/C_functions.c
extern "C" {
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
                          SEXP map_L3_to_L2_down);
}


//' @title Permutation Tests for Directional Gene Sets
//'
//' @description Calculate permutation enrichment scores for directional gene
//'   sets and update vectors needed to calculate NES and p-values.
//'
//' @inheritParams calc_ES_perm
//' @param L3_m integer vector of the unique number of up-regulated genes found
//'   in the directional gene sets, sorted in ascending order.
//' @param L3_w integer vector; the differences between the total number of
//'   genes and the elements of \code{L3_m} (number of genes not up-regulated in
//'   the set).
//' @param L3_m_down integer vector of the unique number of down-regulated genes
//'   found in the directional gene sets, sorted in ascending order.
//' @param L3_w_down integer vector; the differences between the total number of
//'   genes and the elements of \code{L3_m_down} (number of genes not
//'   down-regulated in the set).
//' @param map_L3_to_L2 a 1-based integer vector that maps the unique number of
//'   up-regulated genes to the unique pairs of up- and down-regulated genes.
//' @param map_L3_to_L2_down a 1-based integer vector that maps the unique
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
void calc_ES_perm_dir(IntegerVector& n_same_sign,
                      IntegerVector& n_as_extreme,
                      NumericVector& sum_ES_perm,
                      const Nullable<IntegerVector> seed,
                      const int nperm,
                      int batch_size,
                      const NumericVector& ES,
                      const IntegerVector& ES_end,
                      const NumericVector& y,
                      const NumericVector& r,
                      const int max_set_size,
                      const double sum_ranks,
                      const IntegerVector& L3_m,
                      const IntegerVector& L3_w,
                      const IntegerVector& L3_m_down,
                      const IntegerVector& L3_w_down,
                      IntegerVector& map_L3_to_L2,
                      IntegerVector& map_L3_to_L2_down)
{
  const int n_genes = y.length(); // number of gene-level statistics
  const int n_pairs = map_L3_to_L2.length(); // number of unique pairs
  const int n_sizes_up = L3_m.length();
  const int n_sizes_down = L3_m_down.length();

  // Vectors updated by C_calc_ES_perm_dir (n_as_extreme is also updated)
  IntegerVector n_perm_neg(n_pairs, 0);
  NumericVector sum_perm_pos(n_pairs, 0.0);
  NumericVector sum_perm_neg(n_pairs, 0.0);

  // Convert from 1-based to 0-based indices
  for (int i = 0; i < n_pairs; ++i) {
    --map_L3_to_L2[i];
    --map_L3_to_L2_down[i];
  }

  // The NumericVector versions of L3_w and L3_w_down are needed later to divide
  // quantities when calculating permutation ES anyway, so may as well take
  // their inverses here so that multiplication can be used instead of division,
  // since multiplication takes fewer CPU cycles.
  NumericVector inv_L3_w(n_sizes_up);

  for (int i = 0; i < n_sizes_up; ++i) {
    inv_L3_w[i] = 1.0 / (double)(L3_w[i]);
  }

  NumericVector inv_L3_w_down(n_sizes_down);

  for (int i = 0; i < n_sizes_down; ++i) {
    inv_L3_w_down[i] = 1.0 / (double)(L3_w_down[i]);
  }

  IntegerVector ES_start(n_pairs, 0);

  for (int i = 1; i < n_pairs; ++i) {
    ES_start[i] = ES_end[i - 1];
  }

  // Index of first positive ES by unique combination of the number of
  // up-regulated and down-regulated genes
  IntegerVector ES_pos_idx = get_ES_pos_idx(n_pairs, ES, ES_start, ES_end);

  // Starting positions for segmented sums (up-regulated genes)
  IntegerVector start_vec_up(n_sizes_up, 0);

  for (int i = 1; i < n_sizes_up; ++i) {
    start_vec_up[i] = L3_m[i - 1];
  }

  // Starting positions for segmented sums (down-regulated genes). Avoids
  // collisions with up-regulated genes.
  IntegerVector start_vec_down(n_sizes_down);
  IntegerVector end_vec_down(n_sizes_down);

  start_vec_down[0] = max_set_size - 1;
  end_vec_down[0] = max_set_size - L3_m_down[0] - 1;

  for (int i = 1; i < n_sizes_down; ++i) {
    start_vec_down[i] = end_vec_down[i - 1];
    end_vec_down[i] = max_set_size - L3_m_down[i];
  }

  int n_batches = std::ceil((double)nperm / (double)batch_size);

  // The first batch may be smaller than the rest if nperm is not an integer
  // multiple of n_batches
  const int leftover_perm = nperm % n_batches;

  dqrng::dqset_seed(seed);
  IntegerMatrix perm_idx(max_set_size, batch_size);

  if (leftover_perm > 0) {
    // Each column is a random sample of size max_set_size from integers 0 to
    // n_genes - 1. Used to select pairs of elements from the vectors y and r.
    for (int j = 0; j < leftover_perm; ++j) {
      perm_idx.column(j) = dqrng::dqsample_int(n_genes, max_set_size);
    }

    C_calc_ES_perm_dir(
      n_as_extreme,
      n_perm_neg,
      sum_perm_pos,
      sum_perm_neg,
      ES,
      ES_start,
      ES_end,
      ES_pos_idx,
      y,
      r,
      perm_idx,
      max_set_size,
      leftover_perm, // batch_size
      sum_ranks,
      start_vec_up,
      L3_m,
      inv_L3_w,
      start_vec_down,
      end_vec_down,
      inv_L3_w_down,
      map_L3_to_L2,
      map_L3_to_L2_down
    );

    --n_batches;
  }

  // Loop over remaining permutations (integer multiple of n_batches)
  for (int batch = 0; batch < n_batches; ++batch) {
    for (int j = 0; j < batch_size; ++j) {
      perm_idx.column(j) = dqrng::dqsample_int(n_genes, max_set_size);
    }

    C_calc_ES_perm_dir(
      n_as_extreme,
      n_perm_neg,
      sum_perm_pos,
      sum_perm_neg,
      ES,
      ES_start,
      ES_end,
      ES_pos_idx,
      y,
      r,
      perm_idx,
      max_set_size,
      batch_size,
      sum_ranks,
      start_vec_up,
      L3_m,
      inv_L3_w,
      start_vec_down,
      end_vec_down,
      inv_L3_w_down,
      map_L3_to_L2,
      map_L3_to_L2_down
    );
  }

  // Update n_same_sign and sum_ES_perm. n_as_extreme was updated by
  // C_calc_ES_perm_dir.
  update_results(n_same_sign,
                 sum_ES_perm,
                 n_pairs,
                 nperm,
                 ES_start,
                 ES_end,
                 ES_pos_idx,
                 n_perm_neg,
                 sum_perm_pos,
                 sum_perm_neg);
}
