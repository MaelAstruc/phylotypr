#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp:plugins(openmp)]]

// [[Rcpp::export]]
NumericMatrix calculate_log_probability_cpp3(const NumericMatrix& kmer_genus_count,
                                             const NumericVector& word_specific_priors,
                                             const NumericVector& genus_counts){

  int n_kmers = kmer_genus_count.rows();
  int n_genera = kmer_genus_count.cols();

  NumericMatrix log_probability(n_kmers, n_genera);

  // log((kmer_genus_count + word_specific_priors) / (genus_counts + 1)

  #pragma omp parallel for
  for(int j = 0; j<n_genera; j++){
    log_probability.column(j) = log(
      (kmer_genus_count.column(j) + word_specific_priors)/(genus_counts[j] + 1)
    );
  }
  return log_probability;
}
