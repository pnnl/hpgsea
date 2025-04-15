// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;


//' @title Fast Vector Indexing
//'
//' @param x a numeric or integer vecctor.
//' @param idx an integer vector of indices of \code{x}. May be longer than
//'   \code{x}.
//'
//' @returns The result of \code{x[idx]}.
//'
//' @note These functions are more than twice as fast as their R counterparts.
//'
//' @noRd
//'
// [[Rcpp::export(.Rcpp_indexNumericVector)]]
NumericVector Rcpp_indexNumericVector(NumericVector x, IntegerVector idx) {
   int n = idx.size();

   NumericVector out(n);

   for (int i = 0; i < n; i++) {
      out[i] = x[idx[i] - 1];
   }

   return out;
}

// [[Rcpp::export(.Rcpp_indexIntegerVector)]]
IntegerVector Rcpp_indexIntegerVector(IntegerVector x, IntegerVector idx) {
   int n = idx.size();

   IntegerVector out(n);

   for (int i = 0; i < n; i++) {
      out[i] = x[idx[i] - 1];
   }

   return out;
}


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
//'
// [[Rcpp::export(.Rcpp_matmult_dense)]]
arma::mat Rcpp_matmult_dense(arma::mat& X, arma::mat& Y) {
   arma::mat Z = X * Y;

   return Z;
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
//'
// [[Rcpp::export(.Rcpp_matmult_sparse)]]
arma::mat Rcpp_matmult_sparse(arma::mat& X, arma::sp_mat& Y) {
   arma::mat Z = X * Y;

   return Z;
}


//' @title Core of the .calcES R function
//'
//' @param alpha numeric (\eqn{\geq 0}); the power to which the absolute values
//'   of the entries of \code{X} were raised to construct \code{absStat}. If
//'   \code{alpha=0}, computation time may be significantly reduced, though all
//'   genes in each set will contribute equally.
//' @param Y absolute values of the matrix \code{X} raised to the power of
//'   \code{alpha}. Missing values are then imputed with 0.
//' @param R matrix of ranks of the values in each column of \code{X}.
//'   Missing values in \code{X} are assigned a rank of \code{NA}, which are
//'   then imputed with 0.
//' @param sumRanks integer vector; the sums of the ranks in each sample. Equal
//'   to \code{colSums(rankMat)}.
//' @param A sparse incidence matrix with gene sets as rows and genes as
//'   columns. A value of 1 indicates that the gene is an element of the set,
//'   while a value of 0 indicates otherwise.
//' @param M matrix with gene sets as rows and samples as columns, where each
//'   entry is the number of genes with nonmissing values in each set.
//' @param W matrix with the same dimensions as \code{m} where each entry is
//'   the number of genes with nonmissing values \emph{not} in each set.
//'
//' @returns A matrix of real-valued enrichment scores with gene sets as rows
//'   and samples as columns. May contain missing values if the corresponding
//'   entry of \code{m} is less than 2.
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
//'
// [[Rcpp::export(.Rcpp_calcESCore)]]
arma::mat Rcpp_calcESCore(const double alpha,
                          const int min_size,
                          arma::mat& Y,
                          arma::mat& R,
                          arma::colvec& sumRanks,
                          arma::sp_mat& A,
                          arma::mat& M,
                          arma::mat& W)
{
   arma::mat RA = R * A;

   arma::mat up(M.n_rows, M.n_cols, arma::fill::zeros);

   if (alpha == 0.0) {
      up = RA / M;
   } else {
      // % is Hadamard product, * is dot product, / is Hadamard division
      up = ((R % Y) * A) / (Y * A);
   }

   // Subtract the total sum of ranks in each sample from each row of R * A
   RA.each_col() -= sumRanks;

   arma::mat down = RA / W;

   arma::mat ES = up + down;

   // If the set has fewer than min_size elements with nonmissing values, the ES
   // will be 0.
   arma::uvec indices = arma::find(M < min_size);
   ES(indices).zeros();

   return ES;
}


//' @title Core of the .calcESPerm R function
//'
//' @param alpha non-negative real value.
//' @param Y_perm matrix of absolute values of the input matrix \code{X} (see
//'   \code{\link{fast_ssgsea}}) raised to the power \code{alpha}. The number of
//'   rows is the size of the largest gene set, while the number of columns is
//'   the number of permutations. Each column is a random sample of values from
//'   the i-th row of matrix \code{Y} (see \code{Rcpp_calcESCore}).
//' @param R_perm matrix with the same dimensions as \code{Y_perm} containing
//'   the corresponding ranks of the values of the genes that were selected for
//'   \code{Y_perm}.
//' @param sumRanks_i integer; sum of the ranks of all genes for sample i.
//' @param A_perm dense incidence matrix where the number of rows is the number
//'   of unique gene set sizes and the number of columns is the size of the
//'   largest gene set. Indicates which genes to use to calculate the
//'   permutation ES.
//' @param theta_m_i integer vector of unique gene set sizes.
//' @param theta_w_i integer vector of unique number of genes not in each set.
//'
//' @returns A dense matrix of permutation enrichment scores. Rows correspond to
//'   unique gene set sizes and columns to independent permutations.
//'
//' @details This function sacrifices some precision by using floats instead of
//'   doubles, but the tradeoff is that the runtime is halved. Using floats does
//'   not appreciably impact the NES or p-values in the final results.
//'
//' @author Tyler Sagendorf
//'
//' @references Sanderson, C., & Curtin, R. (2016). Armadillo: A template-based
//'   C++ library for linear algebra. The Journal of Open Source Software, 1(2),
//'   26. \url{https://doi.org/10.21105/joss.00026}
//'
//' @noRd
//'
// [[Rcpp::export(.Rcpp_calcESPermCore)]]
arma::fmat Rcpp_calcESPermCore(const double alpha,
                               arma::fmat& Y_perm,
                               arma::fmat& R_perm,
                               const double sumRanks_i,
                               arma::fmat& A_perm,
                               arma::fvec& theta_m_i,
                               arma::fvec& theta_w_i) {
   // Reciprocals of the size vectors
   arma::fvec m_i_inv = 1.0 / theta_m_i;
   arma::fvec w_i_inv = 1.0 / theta_w_i;

   arma::fmat up_perm(A_perm.n_rows, Y_perm.n_cols, arma::fill::zeros);

   arma::fmat AR_perm = A_perm * R_perm;

   if (alpha == 0.0) {
      // Multiply the diagonal matrix of reciprocals of m_j by
      // A_R_perm. Equivalent to dividing each column of
      // A_R_perm by m_j.
      up_perm = arma::diagmat(m_i_inv) * AR_perm;
   } else {
      // * is dot product, % is Hadamard product, / is Hadamard division
      up_perm = (A_perm * (Y_perm % R_perm)) / (A_perm * Y_perm);
   }

   // Multiply the diagonal matrix of reciprocals of w_i by the matrix of
   // negative sums of ranks of genes not in each set. Equivalent to dividing
   // each column of AR_perm by w_i. down_perm is negative, since sumRanks_i >
   // AR_perm.
   arma::fmat down_perm = arma::diagmat(w_i_inv) * (AR_perm - sumRanks_i);

   arma::fmat ES_perm = up_perm + down_perm;

   return ES_perm;
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
int findFirstPositiveIndex(std::vector<double>& x)
{
   const int SIZE_X = x.size();

   int low = 0;
   int mid;
   int high = SIZE_X - 1;

   while (low <= high) {
      mid = (low + high) / 2; // x is never large enough to cause overflow

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
//' @param x sorted vector of true enrichment scores. Missing values not
//'   allowed.
//' @param y vector of permutation enrichment scores (not necessarily sorted).
//'   All values may be missing.
//'
//' @returns A named list with 3 components, each vectors with length
//'   `length(x)`:
//'
//' \describe{
//'   \item{"n_same_sign_b"}{integer vector; the number of permutation ES in
//'   \code{y} with the same sign as the corresponding ES in \code{x}.}
//'   \item{"n_as_extreme_b"}{integer vector; the number of permutation ES in
//'   \code{y} that were at least as extreme as the corresponding ES in
//'   \code{x}. At most \code{NSameSign.b}.}
//'   \item{"sum_ES_perm_b"}{numeric vector; the absolute value of the sum of
//'   the permutation ES in \code{y} that have the same sign as the
//'   corresponding ES in \code{x}.}
//' }
//'
//' @author Tyler Sagendorf
//'
//' @noRd
//'
// [[Rcpp::export(.Rcpp_extractPermInfo)]]
List Rcpp_extractPermInfo(std::vector<double>& x,
                          std::vector<double>& y)
{
   const int SIZE_X = x.size();
   const int SIZE_Y = y.size();

   // Number of values of y with the same sign as each value of x
   std::vector<int> n_same_sign(SIZE_X, 0);

   // Number of values of y that are at least as extreme as each value of x
   std::vector<int> n_as_extreme(SIZE_X, 0);

   // Vector to store values of sum_y_neg or sum_y_pos for each value of x
   std::vector<double> sum_ES_perm(SIZE_X, 0.0);

   // Number of values of y that are negative
   int n_neg_y = 0;

   // Absolute values of the sums of the negative and positive values of y
   double sum_y_neg = 0.0;
   double sum_y_pos = 0.0;

   // Index of the first positive element of x. If all elements are negative,
   // returns x.size().
   int x_pos_index = findFirstPositiveIndex(x);

   int i = 0; // index for values of x

   for (int j = 0; j < SIZE_Y; j++) {
      if (y[j] < 0.0) { // y[j] is negative
         n_neg_y++;
         sum_y_neg -= y[j];

         i = 0; // start by checking most negative value of x

         // While y is less negative than x and we have not gone past the
         // last negative value of x, move forward to the next most negative
         // x.
         while (i < x_pos_index && y[j] > x[i]) {
            n_as_extreme[i]--;
            i++;
         }

      } else { // y[j] is positive
         sum_y_pos += y[j];

         i = SIZE_X - 1; // start by checking most positive value of x

         // While y is less positive than x and we have not gone past the
         // first positive value of x, move backward to the next most
         // positive x.
         while (i >= x_pos_index && y[j] < x[i]) {
            n_as_extreme[i]--;
            i--;
         }

      }
   }

   // Number of positive values of y
   int n_pos_y = SIZE_Y - n_neg_y;

   // Use the number of negative y and the absolute value of the sum of the
   // negative y for the results of all negative x.
   for (int i = 0; i < x_pos_index; i++) {
      n_as_extreme[i] += n_neg_y;
      n_same_sign[i] = n_neg_y;
      sum_ES_perm[i] = sum_y_neg;
   }

   // Use the number of positive y and the sum of the positive y for the results
   // of all positive x.
   for (int i = x_pos_index; i < SIZE_X; i++) {
      n_as_extreme[i] += n_pos_y;
      n_same_sign[i] = n_pos_y;
      sum_ES_perm[i] = sum_y_pos;
   }

   // The '_b' stands for 'batch', since y might be a single batch of
   // permutation enrichment scores.
   List out = List::create(
      Named("n_same_sign_b") = n_same_sign,
      Named("n_as_extreme_b") = n_as_extreme,
      Named("sum_ES_perm_b") = sum_ES_perm
   );

   return out;
}


//' @title Generate a dense incidence matrix for permutation testing
//'
//' @param end positive integer vector; specifies the number of entries in each
//'   row of the result matrix that will be 1.
//' @param MAX_SET_SIZE integer; the size of the largest set. Determines the
//'   number of rows of the matrix.
//' @param check logical; whether to check each value of \code{end} to ensure it
//'   is greater than 0.
//'
//' @returns An binary 0/1 integer matrix. The number of rows is the length of
//'   \code{end}, and the number of columns is \code{MAX_SET_SIZE}.
//'
//' @author Tyler Sagendorf
//'
//' @noRd
//'
// [[Rcpp::export(.Rcpp_calcAPerm)]]
arma::umat Rcpp_calcAPerm(arma::uvec& end,
                          const int MAX_SET_SIZE,
                          bool check = false)
{
   const int NCOL = end.size();
   // Matrices are in column-major order, so it is better to fill values by
   // column and then take the transpose at the end.
   arma::umat A_perm(MAX_SET_SIZE, NCOL, arma::fill::zeros);

   if (check) { // only needed for directional sets, since they may be empty
      for (int j = 0; j < NCOL; j++) {
         if (end[j] > 0)
            A_perm.col(j).rows(0, end[j] - 1).fill(1);
      }
   } else {
      for (int j = 0; j < NCOL; j++) {
         A_perm.col(j).rows(0, end[j] - 1).fill(1);
      }
   }

   return A_perm.t();
}
