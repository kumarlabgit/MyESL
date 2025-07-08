
#include "gl_logisticr.hpp"
#include "ai_gl_logr.cpp"
#include <iomanip>
#include <sstream>


RRLogisticR::RRLogisticR(const arma::fmat &features, const arma::frowvec &responses, const arma::mat &weights, double *lambda,
                         std::map<std::string, std::string> slep_opts, const bool intercept)
    : lambda(lambda), intercept(intercept) {
    Train(features, responses, weights, slep_opts, intercept);
}

RRLogisticR::RRLogisticR(const arma::fmat &features, const arma::frowvec &responses, const arma::mat &weights, double *lambda,
                         std::map<std::string, std::string> slep_opts, const arma::rowvec &xval_idxs, int xval_id, const bool intercept)
    : lambda(lambda), intercept(intercept) {
    // subset features and responses according to xval_id and xval_idxs
    arma::uvec indices = arma::find(xval_idxs != xval_id);
    Train(features.rows(indices), responses.elem(indices).t(), weights, slep_opts, intercept);
}

void RRLogisticR::writeModelToXMLStream(std::ofstream &XMLFile) {
    int i_level = 0;
    // std::string XMLString = "";
    XMLFile << std::string(i_level * 8, ' ') + "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>" + "\n";
    XMLFile << std::string(i_level * 8, ' ') + "<model>" + "\n";
    i_level++;
    XMLFile << std::string(i_level * 8, ' ') + "<parameters>" + "\n";
    i_level++;
    XMLFile << std::string(i_level * 8, ' ') + "<n_rows>" + std::to_string(this->parameters.n_cols) + "</n_rows>" + "\n";
    XMLFile << std::string(i_level * 8, ' ') + "<n_cols>" + std::to_string(this->parameters.n_rows) + "</n_cols>" + "\n";
    XMLFile << std::string(i_level * 8, ' ') + "<n_elem>" + std::to_string(this->parameters.n_elem) + "</n_elem>" + "\n";
    // for(int i=0; i<this->parameters.n_cols; i++)
    for (int i = 0; i < this->parameters.n_elem; i++) {
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

void RRLogisticR::writeSparseMappedWeightsToStream(std::ofstream &MappedWeightsFile, std::ifstream &FeatureMap) {
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
    // for(int i=0; i<this->parameters.n_cols; i++)
    for (int i = 0; i < this->parameters.n_elem; i++) {
        std::getline(FeatureMap, line);
        if (this->parameters(i) == 0.0) {
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

arma::frowvec RRLogisticR::Train(const arma::fmat &A, const arma::frowvec &responses, const arma::mat &weights, std::map<std::string, std::string> slep_opts,
                              const bool intercept) {

    sll_opts sll_opts_init;

    if (slep_opts.find("maxIter") != slep_opts.end()) {
        sll_opts_init.maxIter = std::stoi(slep_opts["maxIter"]);
    }
    if (slep_opts.find("init") != slep_opts.end()) {
        sll_opts_init.init = std::stoi(slep_opts["init"]);
    }
    if (slep_opts.find("tFlag") != slep_opts.end()) {
        sll_opts_init.tFlag = std::stoi(slep_opts["tFlag"]);
    }
    if (slep_opts.find("nFlag") != slep_opts.end()) {
        sll_opts_init.nFlag = std::stoi(slep_opts["nFlag"]);
    }
    if (slep_opts.find("rFlag") != slep_opts.end()) {
        sll_opts_init.rFlag = std::stoi(slep_opts["rFlag"]);
    }
    if (slep_opts.find("mFlag") != slep_opts.end()) {
        sll_opts_init.mFlag = std::stoi(slep_opts["mFlag"]);
    }
    if (slep_opts.find("lFlag") != slep_opts.end()) {
        sll_opts_init.lFlag = std::stoi(slep_opts["lFlag"]);
    }
    if (slep_opts.find("q") != slep_opts.end()) {
        sll_opts_init.q = std::stod(slep_opts["q"]);
    }
    if (slep_opts.find("rsL2") != slep_opts.end()) {
        sll_opts_init.rsL2 = std::stod(slep_opts["rsL2"]);
    }
    if (slep_opts.find("tol") != slep_opts.end()) {
        sll_opts_init.tol = std::stod(slep_opts["tol"]);
    }

    if (slep_opts.find("nu") != slep_opts.end()) {
		arma::vec nu;
		nu.load(arma::csv_name(slep_opts["nu"]));
		sll_opts_init.nu = nu;
    }
    if (slep_opts.find("mu") != slep_opts.end()) {
		arma::colvec mu;
		mu.load(arma::csv_name(slep_opts["mu"]));
		sll_opts_init.mu = mu;
    }
    if (slep_opts.find("sWeight") != slep_opts.end()) {
		arma::vec sWeight;
		sWeight.load(arma::csv_name(slep_opts["sWeight"]));
		sll_opts_init.sWeight = sWeight;
    }


    arma::uvec temp_ind(weights.col(1).n_elem + 1, arma::fill::zeros);
    temp_ind.subvec(1,weights.col(1).n_elem) = arma::conv_to<arma::uvec>::from(weights.col(1));
    sll_opts_init.ind = temp_ind;

    glLogisticR_result rrModel = glLogisticR(arma::conv_to<arma::mat>::from(A), arma::conv_to<arma::vec>::from(responses), this->lambda[0], sll_opts_init);
    this->parameters = arma::conv_to<arma::fvec>::from(rrModel.x);
    this->nz_gene_count = countNonZeroGenes(parameters, weights);
    return arma::conv_to<arma::frowvec>::from(rrModel.x);

}

/*
    this->intercept = intercept;
    auto trim = [](std::string &s) {
        size_t p = s.find_first_not_of(" \t\r\n");
        s.erase(0, p);

        p = s.find_last_not_of(" \t\r\n");
        if (std::string::npos != p)
            s.erase(p + 1);
    };

    // Set all optional parameters to defaults
    int opts_maxIter = 100;
    int opts_init = 0;
    int opts_tFlag = 5;
    int opts_nFlag = 0;
    int opts_rFlag = 1;
    int opts_mFlag = 0;
    int opts_q = 2;
    double opts_tol = 0.0001;
    int opts_disableEC = 0;
    arma::mat opts_ind = weights;

    // Overwrite default options with those found in slep_opts file.
    if (slep_opts.find("maxIter") != slep_opts.end()) {
        opts_maxIter = std::stoi(slep_opts["maxIter"]);
    }
    int opts_rStartNum = opts_maxIter;
    if (slep_opts.find("init") != slep_opts.end()) {
        opts_init = std::stoi(slep_opts["init"]);
    }
    if (slep_opts.find("tFlag") != slep_opts.end()) {
        opts_tFlag = std::stoi(slep_opts["tFlag"]);
    }
    if (slep_opts.find("nFlag") != slep_opts.end()) {
        opts_nFlag = std::stoi(slep_opts["nFlag"]);
    }
    if (slep_opts.find("rFlag") != slep_opts.end()) {
        opts_rFlag = std::stoi(slep_opts["rFlag"]);
    }
    if (slep_opts.find("mFlag") != slep_opts.end()) {
        opts_mFlag = std::stoi(slep_opts["mFlag"]);
    }
    if (slep_opts.find("tol") != slep_opts.end()) {
        opts_tol = std::stod(slep_opts["tol"]);
    }
    if (slep_opts.find("q") != slep_opts.end()) {
        opts_q = std::stod(slep_opts["q"]);
    }
    if (slep_opts.find("disableEC") != slep_opts.end()) {
        opts_disableEC = std::stoi(slep_opts["disableEC"]);
    }
    std::string line;
    if (slep_opts.find("nu") != slep_opts.end()) {
        std::vector<float> opts_nu;
        std::ifstream nuFile(slep_opts["nu"]);
        if (nuFile.is_open()) {
            while (getline(nuFile, line)) {
                trim(line);
                opts_nu.push_back(std::stod(line));
            }
        }
    }
    if (slep_opts.find("mu") != slep_opts.end()) {
        std::vector<float> opts_mu;
        std::ifstream muFile(slep_opts["mu"]);
        if (muFile.is_open()) {
            while (getline(muFile, line)) {
                trim(line);
                opts_mu.push_back(std::stod(line));
            }
        }
    }
    std::vector<float> opts_sWeight;
    if (slep_opts.find("sWeight") != slep_opts.end()) {
        std::ifstream sWeightFile(slep_opts["sWeight"]);
        if (sWeightFile.is_open()) {
            while (getline(sWeightFile, line)) {
                trim(line);
                opts_sWeight.push_back(std::stod(line));
            }
        }
    }

    // We want to calculate the a_i coefficients of:
    // \sum_{i=0}^n (a_i * x_i^i)
    // In order to get the intercept value, we will add a row of ones.

    // We store the number of rows and columns of the features.
    // Reminder: Armadillo stores the data transposed from how we think of it,
    //           that is, columns are actually rows (see: column major order).
    const size_t nCols = A.n_rows;

    //  arma::mat p = features;
    //  arma::rowvec r = responses;

    // arma::mat& ind = opts_ind;
    int k = opts_ind.n_cols;
    arma::frowvec &ind = arma::frowvec(k + 1, arma::fill::zeros) ind.tail(k) = opts_ind.row(1);
    arma::frowvec &g_weights = opts_ind.row(1);
    arma::fcolvec y = responses.t();
    double *z;
    z = this->Lambda();
    double lambda2_max;
    const size_t m = A.n_rows;
    const size_t n = A.n_cols;
    arma::fcolvec m_ones(m, arma::fill::ones);
    arma::fcolvec m_zeros(m, arma::fill::zeros);
    arma::fcolvec n_zeros(n, arma::fill::zeros);

    double lambda1 = z[0];
    double lambda2 = z[1];

    if (lambda1 < 0 || lambda2 < 0) {
        throw std::invalid_argument("\n z should be nonnegative!\n");
    }

    if (opts_ind.n_cols != 3) {
        throw std::invalid_argument("\n Check opts_ind, expected 3 cols\n");
    }

    if (ind[opts_ind.n_cols + 1] != k) {
        throw std::invalid_argument("\n Check opts_ind, mismatch between last index and feature count\n");
    }

    arma::fcolvec sample_weights(m);
    arma::uvec p_flag = arma::find(y == 1);
    arma::uvec not_p_flag = arma::find(y != 1);
    double m1, m2;
    if (opts_sWeight.size() == 2) {
        std::cout << "Using sample weights of " << opts_sWeight[0] << "(positive) and " << opts_sWeight[1] << "(negative)" << std::endl;
        m1 = p_flag.n_elem * opts_sWeight[0];
        m2 = not_p_flag.n_elem * opts_sWeight[1];
        sample_weights(p_flag).fill(opts_sWeight[0] / (m1 + m2));
        sample_weights(not_p_flag).fill(opts_sWeight[1] / (m1 + m2));

    } else if (opts_sWeight.size() != 0) {
        std::cout << "Invalid sample weights specified, defaulting to unweighted samples." << std::endl;
        sample_weights.fill(1.0 / m);
    } else {
        sample_weights.fill(1.0 / m);
    }

    arma::fcolvec b(m);
    m1 = arma::sum(sample_weights(p_flag)) / arma::sum(sample_weights);
    m2 = 1 - m1;

    double lambda;

    if (opts_rFlag == 0) {
        lambda = z[0];
    } else {
        if (lambda1 < 0 || lambda1 > 1 || lambda2 < 0 || lambda2 > 1) {
            throw std::invalid_argument("\n opts.rFlag=1, so z should be in [0,1]\n");
        }
        b(p_flag) = arma::fcolvec(p_flag.n_elem, arma::fill::ones) * m2;
        b(not_p_flag) = arma::fcolvec(not_p_flag.n_elem, arma::fill::ones) * (-m1);
        b = b % sample_weights;

        arma::fmat ATb = A.t() * b;

        float q_bar;
        if (q <= 1) {
            throw std::invalid_argument("\n q must be greater than 1\n");
        } else {
            q_bar = q / (q - 1);
        }

        norm_ATb = arma::fcolvec(k, arma::fill::ones);
        for (int i = 0; i < k; i++) {
            norm_ATb[i] = norm(ATb.subvec(opts_ind[i], opts_ind[i + 1] - 1), q_bar);
        }

        norm_ATb = norm_ATb / g_weights;

        lambda_max = max(norm_ATb);
        lambda = z[0] * lambda_max;
    }

    arma::fcolvec x(n, arma::fill::zeros); // x.fill(0);
    double c = std::log(m1 / m2);
    arma::fmat Ax = A * x;

    int bFlag = 0;
    double L = 1.0 / m;

    arma::fcolvec weighty = sample_weights % y;

    arma::fcolvec xp = x;
    arma::fcolvec Axp = Ax;
    arma::fcolvec xxp(n, arma::fill::zeros);
    double cp = c;
    double ccp = 0;

    double alphap = 0;
    double alpha = 1;

    double beta, sc, gc, fun_s, fun_x, l_sum, r_sum, tree_norm;
    //  arma::mat As;
    // arma::colvec aa, bb, v, s, prob;
    arma::fcolvec aa, bb, v, s, prob;
    arma::rowvec ValueL(opts_maxIter);
    arma::rowvec funVal(opts_maxIter);
    // arma::mat ind_work(ind.n_rows + 1, ind.n_cols);

    // std::cout << "3..." << std::endl;

    std::cout << "m:" << m << " n:" << n << std::endl;

    for (int iterStep = 0; iterStep < opts_maxIter; iterStep = iterStep + 1) {
        // sgLogisticR.m:293-304
        // std::cout << "4(" << iterStep << ")..." << std::endl;
        beta = (alphap - 1) / alpha;
        s = x + (xxp * beta);
        sc = c + (beta * ccp);
        // std::cout << "sc:" << sc << " c:" << c << " beta:" << beta << " ccp:" << ccp << " m1:" << m1 << " m2:" << m2 << std::endl;
        // std::cout << "4.1..." << std::endl;
        // std::cout << "As_elems:" << As.n_elem << "Ax_elems:" << Ax.n_elem << " Axp_elems:" << Axp.n_elem << " beta:" << beta << std::endl;
        arma::fmat As = Ax + ((Ax - Axp) * beta);
        // std::cout << "4.2..." << std::endl;
        aa = -y % (As + sc);
        // std::cout << "4.3..." << std::endl;
        // sgLogisticR.m:306-316
        bb = arma::max(aa, m_zeros);
        // std::cout << "5..." << std::endl;
        fun_s = arma::as_scalar(sample_weights.t() * (arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb));
        // std::cout << "5.1..." << std::endl;
        prob = m_ones / (m_ones + arma::exp(aa));
        // std::cout << "5.2..." << std::endl;
        b = -weighty % (m_ones - prob);
        // std::cout << "5.3..." << std::endl;
        gc = arma::sum(b);
        // std::cout << "5.4..." << std::endl;
        // sgLogisticR.m:318-329
        // std::cout << "A_rows:" << A.n_rows << " A_cols:" << A.n_cols << " b_elems:" << b.n_elem << std::endl;
        arma::fmat g = A.t() * b;
        // std::cout << "5.5..." << std::endl;
        xp = x;
        Axp = Ax;
        cp = c;
        // std::cout << "5.6..." << std::endl;
        // sgLogisticR.m:331-378
        // std::cout << "6..." << std::endl;
        while (true) {
            v = s - g / L;
            c = sc - gc / L;
            if (q < 1000000) {
                x = eppVector(v, ind, k, n, lambda / L * g_weights, q);
            } else {
                x = eppVector(v, ind, k, n, lambda / L * g_weights, 1000000);
            }

            v = x - s;
            int nonzero_x_count = 0;
            for (int i = 0; i < x.n_rows; i = i + 1) {
                // std::cout << "aa["<< i << "]:" << aa(i) << " bb[" << i << "]:" << bb(i) << " Ax[" << i << "]:" << Ax(i) << std::endl;
                if (x(i) != 0) {
                    nonzero_x_count = nonzero_x_count + 1;
                }
            }

            Ax = A * x;
            // std::cout << "10..." << std::endl;
            // sgLogisticR.m:356-377
            aa = -y % (Ax + c);
            // std::cout << "11..." << std::endl;
            bb = arma::max(aa, m_zeros);
            // std::cout << "aa.n_rows:" << aa.n_rows << " bb.n_rows:" << bb.n_rows << " c:" << c << " sc:" << sc << " gc:" << gc << " L:" << L << std::endl;
            //       for (int i = 0; i <= bb.n_rows; i = i + 1)
            //       {
            //		 std::cout << "aa["<< i << "]:" << aa(i) << " bb[" << i << "]:" << bb(i) << " Ax[" << i << "]:" << Ax(i) << " y[" << i << "]:" << y(i)
            //<< std::endl;
            //	  }
            fun_x = arma::as_scalar(sample_weights.t() * (arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb));
            // std::cout << "12..." << std::endl;
            r_sum = (arma::as_scalar(v.t() * v) + std::pow(c - sc, 2)) / 2;
            l_sum = fun_x - fun_s - arma::as_scalar(v.t() * g) - (c - sc) * gc;
            // std::cout << "13..." << std::endl;
            // std::cout << "r_sum:" << r_sum << " l_sum:" << l_sum << " L:" << L << " fun_x:" << fun_x << std::endl;

            if (r_sum <= std::pow(0.1, 20)) {
                bFlag = 1;
                break;
            }

            //	  if (l_sum <= r_sum * L) // Changed to epsilon comparison on 4/17/2024 to match windows output
            if ((opts_disableEC == 0 && (l_sum < r_sum * L || abs(l_sum - (r_sum * L)) < std::pow(0.1, 12)) || opts_disableEC == 1 && l_sum <= r_sum * L)) {
                break;
            } else {
                L = std::max(static_cast<double>(2 * L), l_sum / r_sum);
            }
        }
        // std::cout << "14..." << std::endl;
        // sgLogisticR.m:382-401
        alphap = alpha;
        alpha = (1 + std::pow(4 * alpha * alpha + 1.0, 0.5)) / 2.0;

        ValueL(iterStep) = L;

        xxp = x - xp;
        ccp = c - cp;
        funVal(iterStep) = fun_x;

        norm_x_k = arma::frowvec(k, arma::fill::zeros);
        for (int i = 0; i < k; i = i + 1) {
            norm_x_k[i] = norm(x.subvec(opts_ind[i], opts_ind[i + 1] - 1), q);
        }

        funVal(iterStep) = fun_x + (lambda * (norm_x_k.t() * g_weights));

        if (bFlag) {
            break;
        }

        switch (opts_tFlag) {
        case 0:
            if (iterStep >= 1) {
                if (std::abs(funVal(iterStep) - funVal(iterStep - 1)) <= opts_tol * funVal(iterStep - 1)) {
                    bFlag = 1;
                }
            }
            break;
        case 5:
            if (iterStep >= opts_maxIter) {
                bFlag = 1;
            }
            break;
        }
        if (bFlag) {
            break;
        }
    }
}
*/

int countNonZeroGenes(const arma::fvec& arr, const arma::mat& ranges) {
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
