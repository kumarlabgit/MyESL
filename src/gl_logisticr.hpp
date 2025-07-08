#include <armadillo>
#include <stdexcept>


class RRLogisticR
{
 public:

  RRLogisticR(const arma::fmat& features,
				const arma::frowvec& responses,
				const arma::mat& weights,
				double* lambda,
				std::map<std::string, std::string> slep_opts,
				const bool intercept = true);

  RRLogisticR(const arma::fmat& features,
				const arma::frowvec& responses,
				const arma::mat& weights,
				double* lambda,
				std::map<std::string, std::string> slep_opts,
				const arma::rowvec& xval_idxs,
				int xval_id,
				const bool intercept = true);


  arma::frowvec Train(const arma::fmat& features,
						const arma::frowvec& responses,
						const arma::mat& weights,
						std::map<std::string, std::string> slep_opts,
						const bool intercept = true);

  void writeModelToXMLStream(std::ofstream& XMLFile);
  void writeSparseMappedWeightsToStream(std::ofstream& MappedWeightsFile, std::ifstream& FeatureMap);

/*
  const arma::fcolvec eppVector(const arma::fcolvec& v_in,
								const arma::mat& ind_mat,
								const int k,
								const int n,
								const arma::mat& rho_mat
								const double p);

  const arma::fcolvec eppVectorR(const arma::fcolvec& t_in,
								const arma::mat& u_mat,
								const arma::fcolvec& v_in,
								const arma::mat& ind_mat,
								const int n,
								const int k);

  float norm(const arma::fcolvec& v_in,
			const float n);
*/


  double* Lambda() { return lambda; }
  int NonZeroGeneCount() { return nz_gene_count; }

 private:
  //Non-zero gene count
  int nz_gene_count = 0;
  /**
   * The calculated B.
   * Initialized and filled by constructor to hold the least squares solution.
   */
  arma::fvec parameters;

  /**
   * The Tikhonov regularization parameter for ridge regression (0 for linear
   * regression).
   */
  double* lambda;
  double lambda1 = lambda[0];

  //! Indicates whether first parameter is intercept.
  bool intercept;
  double intercept_value;
};

int countNonZeroGenes(const arma::fvec& arr, const arma::mat& ranges);
