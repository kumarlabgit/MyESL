
#include "sg_lasso.hpp"
#include <sstream>
#include <iomanip>


SGLasso::SGLasso(const arma::mat& features,
                                   const arma::rowvec& responses,
                                   const arma::mat& weights,
                                   double* lambda,
                                   std::map<std::string, std::string> slep_opts,
                                   const bool intercept) :
    lambda(lambda),
    intercept(intercept)
{
  Train(features, responses, weights, slep_opts, intercept);
}


SGLasso::SGLasso(const arma::mat& features,
                                   const arma::rowvec& responses,
                                   const arma::mat& weights,
                                   double* lambda,
                                   std::map<std::string, std::string> slep_opts,
                                   const arma::rowvec& xval_idxs,
                                   int xval_id,
                                   const bool intercept) :
    lambda(lambda),
    intercept(intercept)
{
  //subset features and responses according to xval_id and xval_idxs
  arma::uvec indices = arma::find(xval_idxs != xval_id);
  Train(features.rows(indices), responses.elem(indices).t(), weights, slep_opts, intercept);
}


void SGLasso::writeModelToXMLStream(std::ofstream& XMLFile)
{
  int i_level = 0;
  //std::string XMLString = "";
  XMLFile << std::string(i_level * 8, ' ') + "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<model>" + "\n";
  i_level++;
  XMLFile << std::string(i_level * 8, ' ') + "<parameters>" + "\n";
  i_level++;
  XMLFile << std::string(i_level * 8, ' ') + "<n_rows>" + std::to_string(this->parameters.n_cols) + "</n_rows>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<n_cols>" + std::to_string(this->parameters.n_rows) + "</n_cols>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<n_elem>" + std::to_string(this->parameters.n_elem) + "</n_elem>" + "\n";
  //for(int i=0; i<this->parameters.n_cols; i++)
  for(int i=0; i<this->parameters.n_elem; i++)
  {
    std::ostringstream streamObj;
    streamObj << std::setprecision(17) << std::scientific << this->parameters(i);
    XMLFile << std::string(i_level * 8, ' ') + "<item>" + streamObj.str() + "</item>" + "\n";
  }
  i_level--;
  XMLFile << std::string(i_level * 8, ' ') + "</parameters>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<lambda1>" + std::to_string(this->lambda[0]) + "</lambda1>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<lambda2>" + std::to_string(this->lambda[1]) + "</lambda2>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<intercept_value>" + std::to_string(this->intercept_value) + "</intercept_value>" + "\n";
  i_level--;
  XMLFile << std::string(i_level * 8, ' ') + "</model>" + "\n";

}

void SGLasso::writeSparseMappedWeightsToStream(std::ofstream& MappedWeightsFile, std::ifstream& FeatureMap)
{
  /*
  int i_level = 0;
  //std::string XMLString = "";
  XMLFile << std::string(i_level * 8, ' ') + "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<model>" + "\n";
  i_level++;
  XMLFile << std::string(i_level * 8, ' ') + "<parameters>" + "\n";
  i_level++;
  XMLFile << std::string(i_level * 8, ' ') + "<n_rows>" + std::to_string(this->parameters.n_cols) + "</n_rows>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<n_cols>" + std::to_string(this->parameters.n_rows) + "</n_cols>" + "\n";
  XMLFile << std::string(i_level * 8, ' ') + "<n_elem>" + std::to_string(this->parameters.n_elem) + "</n_elem>" + "\n";
  */
  std::string line;
  std::getline(FeatureMap, line);
  //for(int i=0; i<this->parameters.n_cols; i++)
  for(int i=0; i<this->parameters.n_elem; i++)
  {
    std::getline(FeatureMap, line);
    if (this->parameters(i) == 0.0)
    {
	  continue;
    }
    std::istringstream iss(line);
    std::string feature_label;
    std::getline(iss, feature_label, '\t');
    std::getline(iss, feature_label, '\t');
    std::ostringstream streamObj;
    streamObj << std::setprecision(17) << std::scientific << this->parameters(i);
    MappedWeightsFile << feature_label + "	" + streamObj.str() + "\n";
  }
  MappedWeightsFile << "Intercept	" + std::to_string(this->intercept_value) + "\n";
  FeatureMap.clear();
  FeatureMap.seekg(0);
}


arma::rowvec& SGLasso::Train(const arma::mat& A,
                               const arma::rowvec& responses,
                               const arma::mat& weights,
                               std::map<std::string, std::string> slep_opts,
                               const bool intercept)
{
  this->intercept = intercept;
  auto trim = [](std::string& s)
  {
     size_t p = s.find_first_not_of(" \t\r\n");
     s.erase(0, p);

     p = s.find_last_not_of(" \t\r\n");
     if (std::string::npos != p)
        s.erase(p+1);
  };

  //Set all optional parameters to defaults
  int opts_maxIter = 100;
  int opts_init = 0;
  int opts_tFlag = 5;
  int opts_nFlag = 0;
  int opts_rFlag = 1;
  int opts_mFlag = 0;
  double opts_tol = 0.0001;
  int opts_disableEC = 0;
  arma::mat opts_ind = weights;

  //Overwrite default options with those found in slep_opts file.
  if ( slep_opts.find("maxIter") != slep_opts.end() ) {
	opts_maxIter = std::stoi(slep_opts["maxIter"]);
  }
  int opts_rStartNum = opts_maxIter;
  if ( slep_opts.find("init") != slep_opts.end() ) {
	opts_init = std::stoi(slep_opts["init"]);
  }
  if ( slep_opts.find("tFlag") != slep_opts.end() ) {
	opts_tFlag = std::stoi(slep_opts["tFlag"]);
  }
  if ( slep_opts.find("nFlag") != slep_opts.end() ) {
	opts_nFlag = std::stoi(slep_opts["nFlag"]);
  }
  if ( slep_opts.find("rFlag") != slep_opts.end() ) {
	opts_rFlag = std::stoi(slep_opts["rFlag"]);
  }
  if ( slep_opts.find("mFlag") != slep_opts.end() ) {
	opts_mFlag = std::stoi(slep_opts["mFlag"]);
  }
  if ( slep_opts.find("tol") != slep_opts.end() ) {
	opts_tol = std::stod(slep_opts["tol"]);
  }
  if ( slep_opts.find("disableEC") != slep_opts.end() ) {
	opts_disableEC = std::stoi(slep_opts["disableEC"]);
  }
  std::string line;
  if ( slep_opts.find("nu") != slep_opts.end() ) {
        std::vector<float> opts_nu;
        std::ifstream nuFile (slep_opts["nu"]);
        if (nuFile.is_open())
        {
          while (getline(nuFile, line))
          {
            trim(line);
            opts_nu.push_back(std::stod(line));
          }
        }
  }
  if ( slep_opts.find("mu") != slep_opts.end() ) {
        std::vector<float> opts_mu;
        std::ifstream muFile (slep_opts["mu"]);
        if (muFile.is_open())
        {
          while (getline(muFile, line))
          {
            trim(line);
            opts_mu.push_back(std::stod(line));
          }
        }
  }
  std::vector<float> opts_sWeight;
  if ( slep_opts.find("sWeight") != slep_opts.end() ) {
        std::ifstream sWeightFile (slep_opts["sWeight"]);
        if (sWeightFile.is_open())
        {
          while (getline(sWeightFile, line))
          {
            trim(line);
            opts_sWeight.push_back(std::stod(line));
          }
        }
  }

  /*
   * We want to calculate the a_i coefficients of:
   * \sum_{i=0}^n (a_i * x_i^i)
   * In order to get the intercept value, we will add a row of ones.
   */

  // We store the number of rows and columns of the features.
  // Reminder: Armadillo stores the data transposed from how we think of it,
  //           that is, columns are actually rows (see: column major order).
  const size_t nCols = A.n_rows;

//  arma::mat p = features;
//  arma::rowvec r = responses;


  arma::mat& ind = opts_ind;
  arma::colvec y = responses.t();
  double* z;
  z = this->Lambda();
  double lambda2_max;
  const size_t m = A.n_rows;
  const size_t n = A.n_cols;
  arma::colvec m_ones(m, arma::fill::ones);
  arma::colvec m_zeros(m, arma::fill::zeros);
  arma::colvec n_zeros(n, arma::fill::zeros);

  double lambda1 = z[0];
  double lambda2 = z[1];

  if (lambda1<0 || lambda2<0)
  {
	throw std::invalid_argument("\n z should be nonnegative!\n");
  }

  if(ind.n_cols != 3)
  {
	throw std::invalid_argument("\n Check opts_ind, expected 3 cols\n");
  }

  //sgLogisticR.m:177-195
  arma::colvec sample_weights(m);
  arma::uvec p_flag = arma::find(y == 1);
  arma::uvec not_p_flag = arma::find(y != 1);
  double m1, m2;
  if (opts_sWeight.size() == 2)
  {
    std::cout << "Using sample weights of " << opts_sWeight[0] << "(positive) and " << opts_sWeight[1] << "(negative)" << std::endl;
    m1 = p_flag.n_elem * opts_sWeight[0];
    m2 = not_p_flag.n_elem * opts_sWeight[1];
    sample_weights(p_flag).fill(opts_sWeight[0] / (m1 + m2));
    sample_weights(not_p_flag).fill(opts_sWeight[1] / (m1 + m2));

  } else if (opts_sWeight.size() != 0) {
    std::cout << "Invalid sample weights specified, defaulting to unweighted samples." << std::endl;
    sample_weights.fill(1.0/m);
  } else {
  sample_weights.fill(1.0/m);
  }

//std::cout << "1..." << std::endl;

  //sgLogisticR.m:200-202
  arma::colvec b(m);
//std::cout << "p_flag.n_elem:" << p_flag.n_elem << std::endl;
  //double m1 = static_cast<double>(p_flag.n_elem) / (double)m;
  m1 = arma::sum(sample_weights(p_flag)) / arma::sum(sample_weights);
  m2 = 1 - m1;

  //sgLogisticR.m:205-241
  double* lambda;

  if (opts_rFlag == 0)
  {
	  lambda = z;
  } else {
	 if (lambda1<0 || lambda1>1 || lambda2<0 || lambda2>1)
	 {
		throw std::invalid_argument("\n opts.rFlag=1, so z should be in [0,1]\n");
	 }
	 b(p_flag) = arma::colvec(p_flag.n_elem, arma::fill::ones) * m2;   b(not_p_flag) = arma::colvec(not_p_flag.n_elem, arma::fill::ones) * (-m1);
	 b = b % sample_weights;

	 arma::mat ATb = A.t() * b;

	 arma::mat temp = arma::abs(ATb);

	 double lambda1_max = arma::as_scalar(arma::max(temp));

	 lambda1 = lambda1 * lambda1_max;

	 temp = arma::max(temp - lambda1, n_zeros);

	 lambda2_max = computeLambda2Max(temp.t(), n, ind, ind.n_rows);

	 lambda2 = lambda2 * lambda2_max;
std::cout << "lambda1_max: " << lambda1_max << " lambda1: " << lambda1 << std::endl;
std::cout << "lambda2_max: " << lambda2_max << " lambda2: " << lambda2 << std::endl;
  }

  //sgLogisticR.m:243-261
  arma::colvec x(n, arma::fill::zeros);   //x.fill(0);
//std::cout << "4..." << std::endl;
  double c = std::log(m1/m2);


//std::cout << "A_cols:" << A.n_cols << " A_rows:" << A.n_rows << " x_cols:" << x.n_cols << " x_rows:" << x.n_rows << std::endl;

  //sgLogisticR.m:264-271
  arma::mat Ax = A * x;

//std::cout << "Ax_cols:" << Ax.n_cols << " Ax_rows:" << Ax.n_rows << std::endl;

//std::cout << "2..." << std::endl;

  //sgLogisticR.m:277-290
  int bFlag = 0;
  double L = 1.0/m;

  arma::colvec weighty = sample_weights % y;

  arma::colvec xp = x;   arma::colvec Axp = Ax;   arma::colvec xxp(n, arma::fill::zeros);
  double cp = c;   double ccp = 0;

  double alphap = 0;   double alpha = 1;


  //sgLogisticR.m:292-442
  double beta, sc, gc, fun_s, fun_x, l_sum, r_sum, tree_norm;
//  arma::mat As;
  arma::colvec aa, bb, prob, s, v;
  arma::rowvec ValueL(opts_maxIter);
  arma::rowvec funVal(opts_maxIter);
  //arma::mat ind_work(ind.n_rows + 1, ind.n_cols);

//std::cout << "3..." << std::endl;

std::cout << "m:" << m << " n:" << n << std::endl;

  for (int iterStep = 0; iterStep < opts_maxIter; iterStep = iterStep + 1)
  {
    //sgLogisticR.m:293-304
//std::cout << "4(" << iterStep << ")..." << std::endl;
    beta = (alphap - 1)/alpha;   s = x + (xxp * beta);   sc = c + (beta * ccp);
//std::cout << "sc:" << sc << " c:" << c << " beta:" << beta << " ccp:" << ccp << " m1:" << m1 << " m2:" << m2 << std::endl;
//std::cout << "4.1..." << std::endl;
//std::cout << "As_elems:" << As.n_elem << "Ax_elems:" << Ax.n_elem << " Axp_elems:" << Axp.n_elem << " beta:" << beta << std::endl;
    arma::mat As = Ax + ((Ax - Axp) * beta);
//std::cout << "4.2..." << std::endl;
    aa = -y % (As + sc);
//std::cout << "4.3..." << std::endl;
    //sgLogisticR.m:306-316
    bb = arma::max(aa,m_zeros);
//std::cout << "5..." << std::endl;
    fun_s = arma::as_scalar(sample_weights.t() * (arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb));
//std::cout << "5.1..." << std::endl;
    prob = m_ones / (m_ones + arma::exp(aa));
//std::cout << "5.2..." << std::endl;
    b = -weighty % (m_ones - prob);
//std::cout << "5.3..." << std::endl;
    gc = arma::sum(b);
//std::cout << "5.4..." << std::endl;
    //sgLogisticR.m:318-329
//std::cout << "A_rows:" << A.n_rows << " A_cols:" << A.n_cols << " b_elems:" << b.n_elem << std::endl;
    arma::mat g = A.t() * b;
//std::cout << "5.5..." << std::endl;
    xp = x;   Axp = Ax;   cp = c;
//std::cout << "5.6..." << std::endl;
    //sgLogisticR.m:331-378
//std::cout << "6..." << std::endl;
    while (true)
    {
      v = s - g/L; c = sc - gc/L;
//std::cout << "7..." << std::endl;
      //sgLogisticR.m:337-338
      arma::mat ind_work = ind;
      ind_work.col(2) = ind_work.col(2) * (lambda2/L);
      arma::rowvec first_row = {-1, -1, lambda1/L};
      ind_work.insert_rows(0, first_row);
      //ind_work(0,0) = -1;
      //ind_work(0,1) = -1;
      //ind_work(0,2) = lambda1/L;
      //ind_work.submat(1,0,ind_work.n_rows,1) = ind.cols(0,1);
      //ind_work(arma::span(1,ind_work.n_rows),2) = ind.col(2) * (lambda2/L);

//std::cout << "8..." << std::endl;
      //sgLogisticR.m:340-342
      x = altra(v, n, ind_work, ind_work.n_rows);
//      arma::colvec temp_x = altra(v, n, ind_work, ind_work.n_rows);
//std::cout << "temp_x_rows:" << temp_x.n_rows << " temp_x_cols:" << temp_x.n_cols << " x_rows:" << x.n_rows << " x_cols:" << x.n_cols << std::endl;
//std::cout << "9..." << std::endl;
      v = x - s;
      int nonzero_x_count = 0;
      for (int i = 0; i < x.n_rows; i = i + 1)
      {
		 //std::cout << "aa["<< i << "]:" << aa(i) << " bb[" << i << "]:" << bb(i) << " Ax[" << i << "]:" << Ax(i) << std::endl;
		 if (x(i) != 0) { nonzero_x_count = nonzero_x_count + 1;}
	  }
//std::cout << "nonzero_x_count:" << nonzero_x_count << std::endl;

      //sgLogisticR.m:345-353
      Ax = A * x;
//std::cout << "10..." << std::endl;
      //sgLogisticR.m:356-377
      aa = -y % (Ax + c);
//std::cout << "11..." << std::endl;
      bb = arma::max(aa, m_zeros);
//std::cout << "aa.n_rows:" << aa.n_rows << " bb.n_rows:" << bb.n_rows << " c:" << c << " sc:" << sc << " gc:" << gc << " L:" << L << std::endl;
//      for (int i = 0; i <= bb.n_rows; i = i + 1)
//      {
//		 std::cout << "aa["<< i << "]:" << aa(i) << " bb[" << i << "]:" << bb(i) << " Ax[" << i << "]:" << Ax(i) << " y[" << i << "]:" << y(i) << std::endl;
//	  }
      fun_x = arma::as_scalar(sample_weights.t() * (arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb));
//std::cout << "12..." << std::endl;
      r_sum = (arma::as_scalar(v.t() * v) + std::pow(c - sc, 2)) / 2;
      l_sum = fun_x - fun_s - arma::as_scalar(v.t() * g) - (c - sc) * gc;
//std::cout << "13..." << std::endl;
//std::cout << "r_sum:" << r_sum << " l_sum:" << l_sum << " L:" << L << " fun_x:" << fun_x << std::endl;

      if (r_sum <= std::pow(0.1, 20))
      {
	     bFlag = 1;
	     break;
	  }

//	  if (l_sum <= r_sum * L) // Changed to epsilon comparison on 4/17/2024 to match windows output
	  if ((opts_disableEC==0 && (l_sum < r_sum * L || abs(l_sum - (r_sum * L)) < std::pow(0.1, 12)) || opts_disableEC==1 && l_sum <= r_sum * L))
	  {
	     break;
	  } else {
	     L = std::max(2*L, l_sum/r_sum);
	  }
    }
//std::cout << "14..." << std::endl;
    //sgLogisticR.m:382-401
    alphap = alpha;   alpha = (1 + std::pow(4 * alpha * alpha + 1.0, 0.5))/2.0;

    ValueL(iterStep) = L;

    xxp = x - xp;   ccp = c - cp;
    funVal(iterStep) = fun_x;

    //ind_work(0,0) = -1;
    //ind_work(0,1) = -1;
    //ind_work(0,2) = lambda1;
    //ind_work.submat(1,0,ind_work.n_rows,1) = ind.cols(0,1);
    //ind_work(arma::span(1,ind_work.n_rows),2) = ind.col(2) * lambda2;
    arma::mat ind_work = ind;
    ind_work.col(2) = ind_work.col(2) * (lambda2);
    arma::rowvec first_row = {-1, -1, lambda1};
    ind_work.insert_rows(0, first_row);

    tree_norm = treeNorm(x.t(), n, ind_work, ind_work.n_rows);

    funVal(iterStep) = fun_x + tree_norm;

    if (bFlag) {break;}

    //sgLogisticR.m:403-435
    switch (opts_tFlag)
    {
	  case 0:
        if (iterStep >=1)
        {
	      if (std::abs(funVal(iterStep) - funVal(iterStep - 1)) <= opts_tol * funVal(iterStep - 1))
	      {
	        bFlag = 1;
	      }
	    }
	    break;
	  case 5:
        if (iterStep >= opts_maxIter)
        {
          bFlag = 1;
        }
        break;
    }
    if (bFlag) {break;}

	//sgLogisticR.m:438-441
	if ((iterStep+1) % opts_rStartNum == 0)
	{
	  alphap = 0;   alpha = 1;   xp = x;   Axp = Ax;   xxp = n_zeros;   L = L/2;
	}

  }

  arma::rowvec x_row = x.as_row();

  parameters = x_row.t();

  std::cout << "Intercept: " << c << std::endl;
  this->intercept_value = c;

//std::cout << "15..." << std::endl;
  this->nz_gene_count = countNonZeroGenes(parameters, weights);

  return x_row;



  // Here we add the row of ones to the features.
  // The intercept is not penalized. Add an "all ones" row to design and set
  // intercept = false to get a penalized intercept.
//  if (intercept)
//  {
//    p.insert_rows(0, arma::ones<arma::mat>(1, nCols));
//  }

//  if (weights.n_elem > 0)
//  {
//    p = p * diagmat(sqrt(weights));
//    r = sqrt(weights) % responses;
//  }

  // Convert to this form:
  // a * (X X^T) = y X^T.
  // Then we'll use Armadillo to solve it.
  // The total runtime of this should be O(d^2 N) + O(d^3) + O(dN).
  // (assuming the SVD is used to solve it)
//  arma::mat cov = p * p.t() +
//      lambda * arma::eye<arma::mat>(p.n_rows, p.n_rows);

//  parameters = arma::solve(cov, p * r.t());
//  return ComputeError(features, responses);

}




/*
void SGLasso::Predict(const arma::mat& points,
    arma::rowvec& predictions) const
{
  if (intercept)
  {
    // We want to be sure we have the correct number of dimensions in the
    // dataset.
    Log::Assert(points.n_rows == parameters.n_rows - 1);
    // Get the predictions, but this ignores the intercept value
    // (parameters[0]).
    predictions = arma::trans(parameters.subvec(1, parameters.n_elem - 1))
        * points;
    // Now add the intercept.
    predictions += parameters(0);
  }
  else
  {
    // We want to be sure we have the correct number of dimensions in
    // the dataset.
    Log::Assert(points.n_rows == parameters.n_rows);
    predictions = arma::trans(parameters) * points;
  }
}

double SGLasso::ComputeError(const arma::mat& features,
                                      const arma::rowvec& responses) const
{
  // Get the number of columns and rows of the dataset.
  const size_t nCols = features.n_cols;
  const size_t nRows = features.n_rows;

  // Calculate the differences between actual responses and predicted responses.
  // We must also add the intercept (parameters(0)) to the predictions.
  arma::rowvec temp;
  if (intercept)
  {
    // Ensure that we have the correct number of dimensions in the dataset.
    if (nRows != parameters.n_rows - 1)
    {
      Log::Fatal << "The test data must have the same number of columns as the "
          "training file." << std::endl;
    }
    temp = responses - (parameters(0) +
        arma::trans(parameters.subvec(1, parameters.n_elem - 1)) * features);
  }
  else
  {
    // Ensure that we have the correct number of dimensions in the dataset.
    if (nRows != parameters.n_rows)
    {
      Log::Fatal << "The test data must have the same number of columns as the "
          "training file." << std::endl;
    }
    temp = responses - arma::trans(parameters) * features;
  }
  const double cost = arma::dot(temp, temp) / nCols;

  return cost;
}
*/

const arma::colvec SGLasso::altra(const arma::colvec& v_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
	double *x;
	x = (double*) malloc(n*sizeof(double));
	const double* v = v_in.memptr();
    int i, j, m;
    double lambda,twoNorm, ratio;
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

//	for(int k=0;k<ind_mat.n_cols*ind_mat.n_rows;k++)
//	{
//	    std::cout << "ind[" << k << "]:" << ind[k] << std::endl;
//	}


//std::cout << "Altra 2..." << std::endl;
    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){

        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }

        lambda=ind[2];

        for(j=0;j<n;j++){
            if (v[j]>lambda)
                x[j]=v[j]-lambda;
            else
                if (v[j]<-lambda)
                    x[j]=v[j]+lambda;
                else
                    x[j]=0;
        }

        i=1;
    }
    else{
        memcpy(x, v, sizeof(double) * n);
        i=0;
    }
//std::cout << "Altra 3..." << std::endl;
    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        lambda=ind[3*i+2];
        if (twoNorm>lambda){
            ratio=(twoNorm-lambda)/twoNorm;

            /*
             * shrinkage this group by ratio
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]*=ratio;
        }
        else{
            /*
             * threshold this group to zero
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]=0;
        }
	}
//std::cout << "Altra 4..." << std::endl;
	arma::colvec x_col(&x[0], n);
	free(x);
//std::cout << "Altra 5..." << std::endl;
	return x_col;
}


const double SGLasso::treeNorm(const arma::rowvec& x_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
	double tree_norm;
	const double* x = x_in.memptr();
    int i, j, m;
    double twoNorm, lambda;
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

//	for(int k=0;k<ind_mat.n_cols*ind_mat.n_rows;k++)
//	{
//	    std::cout << "ind[" << k << "]:" << ind[k] << std::endl;
//	}

    tree_norm=0;

    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){

        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }

        lambda=ind[2];

        for(j=0;j<n;j++){
            tree_norm+=fabs(x[j]);
        }

        tree_norm=tree_norm * lambda;

        i=1;
    }
    else{
        i=0;
    }

    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        lambda=ind[3*i+2];

        tree_norm=tree_norm + lambda*twoNorm;
	}

	return tree_norm;
}


const double SGLasso::computeLambda2Max(const arma::rowvec& x_in,
                            const int n,
                            const arma::mat& ind_mat,
                            const int nodes) const
{
    int i, j, m;
    double lambda,twoNorm;
    const double* x = x_in.memptr();
    double ind[ind_mat.n_cols * ind_mat.n_rows];

    double lambda2_max = 0;

    for(int k=0;k<ind_mat.n_cols;k++)
    {
	    for(int l=0;l<ind_mat.n_rows;l++)
	    {
		    ind[(l*3)+k] = ind_mat(l,k);
		}
	}

    for(i=0;i < nodes; i++){
        /*
         * compute the L2 norm of this group
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];
        twoNorm=sqrt(twoNorm);

        twoNorm=twoNorm/ind[3*i+2];

        if (twoNorm >lambda2_max )
            lambda2_max=twoNorm;
	}

	return lambda2_max;
}


int countNonZeroGenes(const arma::vec& arr, const arma::mat& ranges) {
    auto detectNonZeroInRange = [&arr](int start, int end) -> int {
        for (int i = start; i <= end; ++i) {
            if (arr(i) != 0) {
                return 1;
            }
        }
        return 0;
    };
    int count = 0;

    for (arma::uword i = 0; i < ranges.n_rows; ++i) {
        int start = static_cast<int>(ranges(i, 0))-1;
        int end = static_cast<int>(ranges(i, 1))-1;
        count = count + detectNonZeroInRange(start, end);
    }

    //std::cout << "Number of non-zero genes: " << count << std::endl;
    return count;
}