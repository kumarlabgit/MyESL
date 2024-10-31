
#include <armadillo>
#include <stdexcept>

//namespace mlpack {
//namespace regression /** Regression methods. */ {

/**
 * A simple linear regression algorithm using ordinary least squares.
 * Optionally, this class can perform ridge regression, if the lambda parameter
 * is set to a number greater than zero.
 */
class SGLasso
{
 public:
  /**
   * Creates the model.
   *
   * @param features X, matrix of data points.
   * @param responses y, the measured data for each point in X.
   * @param lambda Regularization constant for ridge regression.
   * @param intercept Whether or not to include an intercept term.
   */
//  SGLasso(const arma::mat& features,
//                   const arma::rowvec& responses,
//                   const double lambda = 0,
//                   const bool intercept = true);

  /**
   * Creates the model with weighted learning.
   *
   * @param features X, matrix of data points.
   * @param responses y, the measured data for each point in X.
   * @param weights Observation weights (for boosting).
   * @param lambda Regularization constant for ridge regression.
   * @param intercept Whether or not to include an intercept term.
   */
  SGLasso(const arma::mat& features,
                   const arma::rowvec& responses,
                   const arma::mat& weights,
                   double* lambda,
                   std::map<std::string, std::string> slep_opts,
                   const bool intercept = true);

  SGLasso(const arma::mat& features,
                   const arma::rowvec& responses,
                   const arma::mat& weights,
                   double* lambda,
                   std::map<std::string, std::string> slep_opts,
                   const arma::rowvec& xval_idxs,
                   int xval_id,
                   const bool intercept = true);

  /**
   * Empty constructor.  This gives a non-working model, so make sure Train() is
   * called (or make sure the model parameters are set) before calling
   * Predict()!
   */
  SGLasso() : lambda(), intercept(true) { }

  /**
   * Train the SGLasso model on the given data. Careful! This will
   * completely ignore and overwrite the existing model. This particular
   * implementation does not have an incremental training algorithm.  To set the
   * regularization parameter lambda, call Lambda() or set a different value in
   * the constructor.
   *
   * @param features X, the matrix of data points to train the model on.
   * @param responses y, the responses to the data points.
   * @param intercept Whether or not to fit an intercept term.
   * @return The least squares error after training.
   */
//  double Train(const arma::mat& features,
//               const arma::rowvec& responses,
//               const bool intercept = true);

  /**
   * Train the SGLasso model on the given data and weights. Careful!
   * This will completely ignore and overwrite the existing model. This
   * particular implementation does not have an incremental training algorithm.
   * To set the regularization parameter lambda, call Lambda() or set a
   * different value in the constructor.
   *
   * @param features X, the matrix of data points to train the model on.
   * @param responses y, the responses to the data points.
   * @param intercept Whether or not to fit an intercept term.
   * @param weights Observation weights (for boosting).
   * @return The least squares error after training.
   */
  arma::rowvec& Train(const arma::mat& features,
               const arma::rowvec& responses,
               const arma::mat& weights,
               std::map<std::string, std::string> slep_opts,
               const bool intercept = true);

  void writeModelToXMLStream(std::ofstream& XMLFile);
  void writeSparseMappedWeightsToStream(std::ofstream& MappedWeightsFile, std::ifstream& FeatureMap);


  /**
   * Calculate y_i for each data point in points.
   *
   * @param points the data points to calculate with.
   * @param predictions y, will contain calculated values on completion.
   */
  //void Predict(const arma::mat& points, arma::rowvec& predictions) const;

  /**
   * Calculate the L2 squared error on the given features and responses using
   * this linear regression model. This calculation returns
   *
   * \f[
   * (1 / n) * \| y - X B \|^2_2
   * \f]
   *
   * where \f$ y \f$ is the responses vector, \f$ X \f$ is the matrix of
   * features, and \f$ B \f$ is the parameters of the trained linear
   * regression model.
   *
   * As this number decreases to 0, the linear regression fit is better.
   *
   * @param points Matrix of features (X).
   * @param responses Transposed vector of responses (y^T).
   */
   /*
  double ComputeError(const arma::mat& points,
                      const arma::rowvec& responses) const;
   */

  const arma::colvec altra(const arma::colvec& v_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const;

  const double treeNorm(const arma::rowvec& x,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const;

  const double computeLambda2Max(const arma::rowvec& x,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const;

  //! Return the parameters (the b vector).
  const arma::vec& Parameters() const { return parameters; }
  //! Modify the parameters (the b vector).
  arma::vec& Parameters() { return parameters; }

  //! Return the Tikhonov regularization parameter for ridge regression.
  double* Lambda() { return lambda; }
  //! Modify the Tikhonov regularization parameter for ridge regression.
  //double& Lambda() { return lambda; }

  //! Return whether or not an intercept term is used in the model.
  bool Intercept() const { return intercept; }

  int NonZeroGeneCount() { return nz_gene_count; }

  /**
   * Serialize the model.
   */

  //template<typename Archive>
  //void serialize(Archive& ar, const unsigned int /* version */)
  /*
  {
    ar & BOOST_SERIALIZATION_NVP(parameters);
    ar & BOOST_SERIALIZATION_NVP(lambda1);
    ar & BOOST_SERIALIZATION_NVP(intercept_value);
  }
  */


 private:
  //Non-zero gene count
  int nz_gene_count = 0;
  /**
   * The calculated B.
   * Initialized and filled by constructor to hold the least squares solution.
   */
  arma::vec parameters;

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

//} // namespace regression
//} // namespace mlpack

//#endif // MLPACK_METHODS_SG_LASSO_HPP

int countNonZeroGenes(const arma::vec& arr, const arma::mat& ranges);
