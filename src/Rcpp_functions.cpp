// [[Rcpp::depends(RcppArmadillo, dqrng)]]
#include <RcppArmadillo.h>
#include <dqrng.h>

using namespace Rcpp;


//' @title Dense Matrix Multiplication
//'
//' @description Multiplication of two dense matrices.
//'
//' @param X,Y dense matrices.
//'
//' @returns The product of \code{X} and \code{Y}: a dense matrix.
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//' @noRd

// [[Rcpp::export(.Cpp_matmult_dense)]]
arma::mat matmult_dense(const arma::mat& X,
                        const arma::mat& Y)
{
  return X * Y;
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
//' @param alpha numeric (\eqn{\geq 0}); the power to which the absolute values
//'   of the entries of \code{X} were raised to construct \code{Y}. If
//'   \code{alpha=0}, computation time may be significantly reduced, though all
//'   genes in each set will contribute equally.
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
arma::mat calcES(const double alpha,
                 const int min_size,
                 const arma::mat& Y,
                 const arma::mat& R,
                 const arma::colvec& sumRanks,
                 const arma::sp_mat& A,
                 const arma::mat& M,
                 const arma::mat& W)
{
  arma::mat RA = R * A;

  arma::mat ES = RA;

  if (alpha == 0.0) {
    ES /= M;
  } else {
    // % is Hadamard product, * is dot product, / is Hadamard division
    ES = ((R % Y) * A) / (Y * A);
  }

  // Subtract the total sum of ranks in each sample from the corresponding row
  // of R * A
  RA.each_col() -= sumRanks;

  ES += RA / W;

  // If the set has fewer than min_size elements with nonmissing values, the ES
  // will be 0. Only applies to directional gene sets.
  arma::uvec indices = arma::find(M < min_size);
  ES(indices).zeros();

  return ES;
}


//' @title Create a Matrix of Random Indices
//'
//' @description Create a matrix where each column is an independent random
//'   sample of indices, taken without replacement from 0 to \code{n_genes - 1}.
//'
//' @param n_genes one more than the maximum index.
//' @param seeds integer vector of random seeds. The length will determine the
//'   number of columns in the output matrix.
//' @param max_set_size the size of the largest gene set. This will determine
//'   the number of rows in the output matrix.
//'
//' @returns A matrix.
//'
//' @seealso
//' https://cran.r-project.org/web/packages/dqrng/vignettes/cpp-api.html
//'
//' @noRd
arma::umat create_perm_indices(const size_t n_genes,
                               const IntegerVector& seeds,
                               const size_t max_set_size)
{
  const size_t batch_size = seeds.size();

  arma::umat perm_indices(max_set_size, batch_size, arma::fill::zeros);

  for (size_t j = 0; j < batch_size; ++j) {
    dqrng::dqset_seed(IntegerVector::create(seeds[j]));

    // Sample a total of max_set_size indices taken without replacement from 0
    // to n_genes - 1
    const IntegerVector sample = dqrng::dqsample_int(n_genes, max_set_size);

    for (size_t i = 0; i < max_set_size; ++i) {
      perm_indices.at(i, j) = sample[i];
    }
  }

  return perm_indices;
}


//' @title Multiplication of an unseen binary matrix and real-valued matrix
//'
//' @description Compute the product of an unseen binary matrix and real-valued
//'   matrix when the binary matrix is in a specific format. See details.
//'
//' @param set_sizes integer vector of unique set sizes, sorted in ascending
//'   order.
//' @param Y real-valued matrix. The number of rows is equal to
//'   \code{max(set_sizes)}.
//'
//' @details The set_sizes vector acts as a stand-in for the unseen binary
//'   matrix, which has the following characteristics:
//'
//'   1. The first \code{set_sizes[i]} elements of the i-th row are all 1 and
//'   the remaining entries of that row are 0.
//'
//'   2. If the i-th row has n ones, then the i + 1 row will have at least (n +
//'   1) ones.
//'
//'   3. All elements of the last row are 1.
//'
//'   This function has time complexity \code{O(max(set_sizes) * Y.n_cols)}.
//'
//' @author Tyler Sagendorf
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//' @noRd
arma::mat binary_matmult(const arma::uvec& set_sizes,
                         const arma::mat& Y)
{
  const arma::uword batch_size = Y.n_cols;
  const arma::uword N_UNIQUE_SIZES = set_sizes.size();

  arma::mat out(N_UNIQUE_SIZES, batch_size, arma::fill::zeros);

  // Define start positions to sum a window of Y column values. End positions
  // are the set sizes.
  arma::uvec start(N_UNIQUE_SIZES, arma::fill::zeros);

  if (N_UNIQUE_SIZES != 1) {
    start.subvec(1, N_UNIQUE_SIZES - 1) =
      set_sizes.subvec(0, N_UNIQUE_SIZES - 2);
  }

  // Loop over permutations
  for (arma::uword j = 0; j < batch_size; ++j) {
    double sum_ij = 0.0; // cumulative sum
    const arma::vec col_j = Y.unsafe_col(j);

    // Loop over unique set sizes
    for (arma::uword i = 0; i < N_UNIQUE_SIZES; ++i) {
      // Sum the values for genes in the set
      for (arma::uword k = start.at(i); k < set_sizes.at(i); ++k) {
        sum_ij += col_j.at(k);
      }

      out.at(i, j) = sum_ij;
    }
  }

  return out;
}


arma::mat binary_matmult_up(const arma::mat& Y,
                            const arma::uvec& set_sizes_up)
{
  const arma::uword batch_size = Y.n_cols;
  const arma::uword N_UNIQUE_SIZES = set_sizes_up.size();

  arma::mat out(N_UNIQUE_SIZES, batch_size, arma::fill::zeros);

  double sum_ij = 0.0;

  // Loop over permutations
  for (arma::uword j = 0; j < batch_size; ++j) {
    const arma::vec col_j = Y.unsafe_col(j);

    // Loop over unique set sizes
    for (arma::uword i = 0; i < N_UNIQUE_SIZES; ++i) {
      sum_ij = 0.0;

      // Sum the values for genes in the set. This also handles empty sets
      for (arma::uword k = 0; k < set_sizes_up.at(i); ++k) {
        sum_ij += col_j.at(k);
      }

      out.at(i, j) = sum_ij;
    }
  }

  return out;
}


arma::mat binary_matmult_down(const arma::mat& Y,
                              const arma::uvec& set_sizes_down)
{
  const arma::uword batch_size = Y.n_cols;
  const arma::uword N_UNIQUE_SIZES = set_sizes_down.size();
  const arma::uword MAX_SET_SIZE = Y.n_rows;

  arma::mat out(N_UNIQUE_SIZES, batch_size, arma::fill::zeros);

  double sum_ij = 0.0;

  // Loop over permutations
  for (arma::uword j = 0; j < batch_size; ++j) {
    const arma::vec col_j = Y.unsafe_col(j);

    // Loop over unique set sizes
    for (arma::uword i = 0; i < N_UNIQUE_SIZES; ++i) {
      sum_ij = 0.0;

      // Sum the values for genes in the set. This also handles empty sets
      for (arma::uword k = MAX_SET_SIZE - set_sizes_down.at(i);
           k < MAX_SET_SIZE;
           ++k) {
        sum_ij += col_j.at(k);
      }

      out.at(i, j) = sum_ij;
    }
  }

  return out;
}


//' @title Calculate Permutation Enrichment Scores
//'
//' @param alpha non-negative real value.
//' @param y_i numeric vector of absolute values raised to the power of
//'   \code{alpha}.
//' @param r_i corresponding vector of ranks of the genes in \code{y_i}.
//' @param seeds integer vector of seeds used for random sampling.
//' @param max_set_size integer; the size of the largest gene set that will be
//'   tested. Determines the size of each random sample.
//' @param sumRanks_i numeric; sum of the ranks of all genes.
//' @param theta_m_i integer vector of unique gene set sizes. These will be
//'   unique and sorted in ascending order.
//' @param theta_w_i integer vector of the unique number of genes not in each
//'   set. Equal to the total number of genes with nonmissing values minus the
//'   elements of \code{theta_m_i}.
//'
//' @returns A dense matrix of permutation enrichment scores. Rows correspond to
//'   unique gene set sizes and columns to independent permutations.
//'
//' @details This function has time complexity \code{O(max_set_size *
//'   seeds.size())}.
//'
//' @author Tyler Sagendorf
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//' @noRd

// [[Rcpp::export(.Cpp_calcESPerm)]]
arma::mat calcESPerm(const double alpha,
                     const arma::vec& y_i,
                     const arma::vec& r_i,
                     const IntegerVector& seeds,
                     const size_t max_set_size,
                     const double sumRanks_i,
                     const arma::uvec& theta_m_i,
                     const arma::vec& theta_w_i)
{
  const size_t n_genes = r_i.size(); // max index for sampling
  const arma::uword batch_size = seeds.size();

  // Each column is an independent random sample of indices
  arma::umat perm_indices = create_perm_indices(n_genes, seeds, max_set_size);

  arma::mat R_perm(max_set_size, batch_size, arma::fill::zeros);

  // Each column of R_perm is an independent random sample of the gene ranks
  for (arma::uword j = 0; j < batch_size; ++j) {
    R_perm.unsafe_col(j) = r_i.elem(perm_indices.unsafe_col(j));
  }

  arma::mat AR_perm = binary_matmult(theta_m_i, R_perm);

  arma::mat ES_perm = AR_perm;

  if (alpha == 0.0) {
    // Divide each column of AR_perm by the unique set sizes.
    ES_perm.each_col() /= arma::conv_to<arma::vec>::from(theta_m_i);
  } else {
    arma::mat Y_perm = R_perm;

    // Each column of Y_perm is an independent random sample of the absolute
    // values of the genes raised to the power of alpha. The ranks in R_perm are
    // the corresponding ranks of the original gene-level values.
    for (arma::uword j = 0; j < batch_size; ++j) {
      Y_perm.unsafe_col(j) = y_i.elem(perm_indices.unsafe_col(j));
    }

    // % is Hadamard product, / is Hadamard division
    ES_perm = binary_matmult(theta_m_i, Y_perm % R_perm) /
      binary_matmult(theta_m_i, Y_perm);
  }

  AR_perm -= sumRanks_i;
  AR_perm.each_col() /= theta_w_i;

  ES_perm += AR_perm;

  return ES_perm;
}


//' @title Calculate Permutation ES for Directional Gene Sets
//'
//' @inheritParams .Cpp_calcESPerm
//'
//' @param theta_m_i integer vector of the number of up-regulated genes in each
//'   set. Each \code{(theta_m_i, theta_m_d_i)} pair is unique.
//' @param theta_w_i integer vector of the number of genes that are not
//'   up-regulated in the set. Will count genes that are either not in the set
//'   or down-regulated in the set.
//' @param theta_m_d_i integer vector of the number of down-regulated genes in
//'   each set. Each \code{(theta_m_i, theta_m_d_i)} pair is unique.
//' @param theta_w_d_i integer vector of the number of genes that are not
//'   down-regulated in the set. Will count genes that are either not in the set
//'   or up-regulated in the set.
//' @param min_size integer; minimum set size. If the number of "up" or "down"
//'   genes in the set is less than \code{min_size}, the directional ES will be
//'   set to 0.
//'
//' @returns A matrix of permutation ES for directional gene sets.
//'
//' @author Tyler Sagendorf
//'
//' @noRd

// [[Rcpp::export(.Cpp_calcESPerm_dir)]]
arma::mat calcESPerm_dir(const double alpha,
                         const arma::vec& y_i,
                         const arma::vec& r_i,
                         const IntegerVector& seeds,
                         const size_t max_set_size,
                         const double sumRanks_i,
                         const arma::uvec& theta_m_i,
                         const arma::vec& theta_w_i,
                         const arma::uvec& theta_m_d_i,
                         const arma::vec& theta_w_d_i,
                         const arma::uword& min_size)
{
  const size_t n_genes = r_i.size(); // max index for sampling
  const arma::uword batch_size = seeds.size();

  // Each column is an independent random sample of indices
  arma::umat perm_indices = create_perm_indices(n_genes, seeds, max_set_size);

  arma::mat R_perm(max_set_size, batch_size, arma::fill::zeros);

  for (arma::uword j = 0; j < batch_size; ++j) {
    R_perm.unsafe_col(j) = r_i.elem(perm_indices.unsafe_col(j));
  }

  arma::mat AR_perm_up = binary_matmult_up(R_perm, theta_m_i);
  arma::mat AR_perm_down = binary_matmult_down(R_perm, theta_m_d_i);

  arma::mat ES_perm_up = AR_perm_up;
  arma::mat ES_perm_down = AR_perm_down;

  if (alpha == 0.0) {
    ES_perm_up.each_col() /= arma::conv_to<arma::vec>::from(theta_m_i);

    ES_perm_down.each_col() /= arma::conv_to<arma::vec>::from(theta_m_d_i);
  } else {
    arma::mat Y_perm = R_perm;

    for (arma::uword j = 0; j < batch_size; ++j) {
      Y_perm.unsafe_col(j) = y_i.elem(perm_indices.unsafe_col(j));
    }

    ES_perm_up = binary_matmult_up(Y_perm % R_perm, theta_m_i) /
        binary_matmult_up(Y_perm, theta_m_i);

    ES_perm_down = binary_matmult_down(Y_perm % R_perm, theta_m_d_i) /
        binary_matmult_down(Y_perm, theta_m_d_i);
  }

  // Average rank of genes, excluding those that are "up" and in the sets
  AR_perm_up -= sumRanks_i;
  AR_perm_up.each_col() /= theta_w_i;
  ES_perm_up += AR_perm_up;

  // Average rank of genes, excluding those that are "down" and in the sets
  AR_perm_down -= sumRanks_i;
  AR_perm_down.each_col() /= theta_w_d_i;
  ES_perm_down += AR_perm_down;

  // If a set is too small or empty, replace all permutation ES in that row with
  // zero.
  const arma::uvec idx_up = arma::find(theta_m_i < min_size);
  const arma::uvec idx_down = arma::find(theta_m_d_i < min_size);

  ES_perm_up.rows(idx_up).zeros();
  ES_perm_down.rows(idx_down).zeros();

  return ES_perm_up - ES_perm_down;
}


//' @title Find Index of First Positive Value in Sorted Vector
//'
//' @description Given a sorted vector of real-valued numbers, find the index of
//'   the first positive value using a binary search.
//'
//' @param x a sorted real-valued vector.
//'
//' @returns The 0-based index (integer) of the first positive value, or the
//'   size of \code{x} if no positive values were found.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
int findFirstPositiveIndex(const std::vector<double>& x)
{
  const int SIZE_X = x.size();

  int low = 0;
  int mid;
  int high = SIZE_X - 1;

  while (low <= high) {
    mid = (low + high) >> 1; // x is never large enough to cause overflow

    if (x[mid] >= 0.0) {
      // check if the value of x is the first positive
      if (mid == 0 || x[mid - 1] < 0.0) {
        return mid;
      } else {
        high = mid - 1;
      }
    } else {
        low = mid + 1;
    }
  }

  return SIZE_X; // no positive values found
}


//' @title Extract Information About Permutation Enrichment Scores
//'
//' @param x sorted vector of true enrichment scores.
//' @param y vector of permutation enrichment scores (not sorted).
//'
//' @returns A \code{data.table} with 3 columns where the number of rows is
//'   equal to the length of \code{x}.
//'
//' \describe{
//'   \item{"n_same_sign_b"}{integer vector; the number of permutation ES in
//'   \code{y} with the same sign as the corresponding ES in \code{x}.}
//'
//'   \item{"n_as_extreme_b"}{integer vector; the number of permutation ES in
//'   \code{y} that were at least as extreme as the corresponding ES in
//'   \code{x}. At most \code{n_same_sign_b}.}
//'
//'   \item{"sum_ES_perm_b"}{numeric vector; the absolute value of the sum of
//'   the permutation ES in \code{y} that have the same sign as the
//'   corresponding ES in \code{x}.}
//' }
//'
//' @author Tyler Sagendorf
//'
//' @noRd
void extractPermInfo_internal(List ls,
                              const std::vector<double>& x,
                              const std::vector<double>& y)
{
  const int SIZE_X = x.size();
  const int SIZE_Y = y.size();

  // Number of values of y that are negative
  int n_neg_y = 0;

  // Absolute values of the sums of the negative and positive values of y
  double sum_y_neg = 0.0;
  double sum_y_pos = 0.0;

  // Number of values of y with the same sign as each value of x
  std::vector<int> n_same_sign = std::vector<int>(ls["n_same_sign"]);

  // Number of values of y that are at least as extreme as each value of x
  std::vector<int> n_as_extreme = std::vector<int>(ls["n_as_extreme"]);

  // Vector to store values of sum_y_neg or sum_y_pos for each value of x
  std::vector<double> sum_ES_perm = std::vector<double>(ls["sum_ES_perm"]);

  // Index of the first positive element of x. If all elements are negative,
  // returns x.size().
  const int x_pos_index = findFirstPositiveIndex(x);

  for (int j = 0; j < SIZE_Y; ++j) {
    const double y_j = y[j]; // cache to avoid repeated indexing

    if (y_j < 0.0) {
      ++n_neg_y;
      sum_y_neg -= y_j;

      for (int i = 0; i < x_pos_index; ++i) {
        n_as_extreme[i] += y_j <= x[i];
      }

    } else {
      sum_y_pos += y_j;

      for (int i = x_pos_index; i < SIZE_X; ++i) {
        n_as_extreme[i] += y_j >= x[i];
      }

    }
  }

  // Number of positive values of y
  const int n_pos_y = SIZE_Y - n_neg_y;

  // Use the number of negative y and the absolute value of the sum of the
  // negative y for the results of all negative x.
  for (int i = 0; i < x_pos_index; ++i) {
    n_same_sign[i] += n_neg_y;
    sum_ES_perm[i] += sum_y_neg;
  }

  // Use the number of positive y and the sum of the positive y for the results
  // of all positive x.
  for (int i = x_pos_index; i < SIZE_X; ++i) {
    n_same_sign[i] += n_pos_y;
    sum_ES_perm[i] += sum_y_pos;
  }

  // Update list vectors. ls is modified by reference
  ls["n_same_sign"] = IntegerVector(n_same_sign.begin(), n_same_sign.end());
  ls["n_as_extreme"] = IntegerVector(n_as_extreme.begin(), n_as_extreme.end());
  ls["sum_ES_perm"] = NumericVector(sum_ES_perm.begin(), sum_ES_perm.end());
}


//' @title Extract Information from a Permutation Enrichment Score Matrix
//'
//' @description Extract information from a matrix of permutation enrichment
//'   scores run as a single batch.
//'
//' @param ES_ls list of sorted true enrichment scores grouped by gene set size.
//' @param ES_perm matrix of permutation ES. The number of rows is equal to the
//'   length of \code{ES_ls}, while the number of columns is at most the total
//'   number of permutations: more likely, it is a fraction of the total number
//'   of permutations. See the \code{batch_size} parameter of
//'   \code{\link{fast_ssgsea}} for more details.
//'
//' @returns A list of \code{data.table} objects, each with 3 columns:
//'
//' \describe{
//'   \item{"n_same_sign_b"}{integer; the number of permutation ES in each
//'   row of \code{ES_perm} with the same sign as the corresponding ES in
//'   \code{ES}.}
//'   \item{"n_as_extreme_b"}{integer; the number of permutation ES in
//'   each row of \code{ES_perm} that were at least as extreme as the
//'   corresponding ES in \code{ES}. At most \code{"n_same_sign_b"}.}
//'   \item{"sum_ES_perm_b"}{integer; the sum of the absolute values of the
//'   permutation ES that have the same sign as the corresponding ES in
//'   \code{ES}.}
//' }
//'
//' @author Tyler Sagendorf
//'
//' @noRd

// [[Rcpp::export(.Cpp_extractPermInfo)]]
void extractPermInfo(List ls,
                     const List ES_ls,
                     const NumericMatrix& ES_perm)
{
  const int N_UNIQUE_SIZES = ES_ls.size();

  NumericVector temp(ES_perm.ncol(), 0.0); // the size doesn't change

  for (int i = 0; i < N_UNIQUE_SIZES; ++i) {
    // Extract i-th row of ES_perm (is there a better way?)
    temp = ES_perm(i, _);
    const std::vector<double> ES_perm_i = as<std::vector<double>>(temp);

    extractPermInfo_internal(ls[i], ES_ls[i], ES_perm_i);
  }
}
