
#include <armadillo>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdexcept>
#include "eppVector.c"
#include "eppVectorR.c"
#include "initFactor.cpp"
#include "sll_opts.cpp"

// Forward declarations - these functions should be imported from reference files

/*
struct sll_opts {
    // Starting point
    std::optional<arma::vec> x0;           // Starting point of x
    std::optional<double> c0;              // Starting point for intercept c
    int init = 0;                          // Initialization method (default: 0)

    // Termination
    int maxIter = 10000;                   // Maximum number of iterations (default: 1e4)
    double tol = 1e-3;                     // Tolerance parameter (default: 1e-4, but code uses 1e-3)
    int tFlag = 0;                         // Flag for termination (default: 0)

    // Normalization
    int nFlag = 0;                         // Flag for implicit normalization of A (default: 0)
    std::optional<arma::rowvec> mu;        // Row vector to be subtracted from each sample
    std::optional<arma::vec> nu;           // Weight vector for normalization

    // Regularization
    int rFlag = 0;                         // Flag for regularization (default: 0)
    double rsL2 = 0.0;                     // Regularization parameter for squared L2 norm (default: 0)

    // Method & Line Search
    int lFlag = 0;                         // Line search flag (default: 0)
    int mFlag = 0;                         // Method flag (default: 0)

    // Group & Others
    std::optional<arma::uvec> ind;         // Indices for k groups
    double q = 2.0;                        // Value of q in L1/Lq regularization (default: 2)
    std::optional<arma::vec> sWeight;      // Sample weight for Logistic Loss
    std::optional<arma::vec> gWeight;      // Weight for different groups
    std::string fName;                     // Name of the function
};
*/


// Note: Function signature mismatch - eppVector expects different parameter types than called with
// arma::vec eppVector(arma::vec& v, const arma::vec& ind, int k, int n, const arma::vec& lambda_over_L_times_gWeight, double q);

// Note: Function signature mismatch - eppVectorR expects different parameter types than called with
// std::pair<arma::vec, arma::vec> eppVectorR(const arma::vec& u, const arma::vec& v, const arma::vec& ind, int n, int k);

//sll_opts sll_opts_default(const sll_opts& input_opts);

struct glLogisticR_result {
    arma::vec x;
    double c;
    arma::vec funVal;
    arma::vec ValueL;
};

glLogisticR_result glLogisticR(const arma::mat& A, const arma::vec& y, double z, const sll_opts& opts) {
    //
    // Function glLogisticR
    //      Logistic Loss with the (group) L1/Lq-norm Regularization.
    //
    //      The features are grouped into k groups.
    //
    // Problem
    //
    //  min  f(x,c) = - sum_i weight_i * log (p_i) + z * gWeight_j sum_j ||x^j||_q
    //
    //  a_i denotes a training sample,
    //      and a_i' corresponds to the i-th row of the data matrix A
    //
    //  y_i (either 1 or -1) is the response
    //
    //  p_i= 1/ (1+ exp(-y_i (x' * a_i + c) ) ) denotes the probability
    //
    //  weight_i denotes the weight for the i-th sample
    //
    //  x is grouped into k groups according to opts.ind.
    //      The indices of x_j in x is (ind(j)+1):ind(j+1).
    //
    // Input parameters:
    //
    //  A-         Matrix of size m x n
    //                A can be a dense matrix
    //                         a sparse matrix
    //                         or a DCT matrix
    //  y -        Response vector (of size mx1)
    //  z -        L1/Lq norm regularization parameter (z >=0)
    //  opts-      Optional inputs (default value: opts=[])
    //             !!For glLogisticR, we require that opts.ind is specified.!!
    //
    // Output parameters:
    //  x-         The obtained weight of size n x 1
    //  c-         The obtained intercept (scalar)
    //  funVal-    Function value during iterations
    //
    // Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
    //
    //
    // Related functions:
    //
    //  sll_opts, initFactor, pathSolutionLeast
    //  glLeastR, eppVectorR
    //

    // Verify and initialize the parameters
    int m = A.n_rows;
    int n = A.n_cols;

    if (y.n_elem != m) {
        throw std::runtime_error("Check the length of y!");
    }

    if (z <= 0) {
        throw std::runtime_error("z should be positive!");
    }

    sll_opts processed_opts = default_sll_opts(opts); // run sll_opts to set default values (flags)

    // Detailed initialization
    // Normalization

    // Please refer to the function 'sll_opts'
    //                 for the definitions of mu, nu and nFlag
    //
    // If .nFlag =1, the input matrix A is normalized to
    //                     A= ( A- repmat(mu, m,1) ) * diag(nu)^{-1}
    //
    // If .nFlag =2, the input matrix A is normalized to
    //                     A= diag(nu)^{-1} * ( A- repmat(mu, m,1) )

    arma::vec mu, nu;
    arma::uvec ind_zero;

    if (processed_opts.nFlag != 0) {
        if (processed_opts.mu.has_value()) {
            mu = processed_opts.mu.value();
            if (mu.n_elem != n) {
                throw std::runtime_error("Check the input .mu");
            }
        } else {
            mu = arma::mean(A, 0).t();
        }

        if (processed_opts.nFlag == 1) {
            if (processed_opts.nu.has_value()) {
                nu = processed_opts.nu.value();
                if (nu.n_elem != n) {
                    throw std::runtime_error("Check the input .nu!");
                }
            } else {
                nu = arma::sqrt(arma::sum(arma::square(A), 0) / m).t();
            }
        } else { // .nFlag=2
            if (processed_opts.nu.has_value()) {
                nu = processed_opts.nu.value();
                if (nu.n_elem != m) {
                    throw std::runtime_error("Check the input .nu!");
                }
            } else {
                nu = arma::sqrt(arma::sum(arma::square(A), 1) / n);
            }
        }

        ind_zero = arma::find(arma::abs(nu) <= 1e-10);
        nu.elem(ind_zero).ones();
        // If some values in nu is typically small, it might be that,
        // the entries in a given row or column in A are all close to zero.
        // For numerical stability, we set the corresponding value to 1.
    }

/*
    if (!A.is_sparse() && processed_opts.nFlag != 0) {
        std::cout << "\n -----------------------------------------------------" << std::endl;
        std::cout << "\n The data is not sparse or not stored in sparse format" << std::endl;
        std::cout << "\n The code still works." << std::endl;
        std::cout << "\n But we suggest you to normalize the data directly," << std::endl;
        std::cout << "\n for achieving better efficiency." << std::endl;
        std::cout << "\n -----------------------------------------------------" << std::endl;
    }
*/

    // Group & Others

    // Initialize ind and q
    arma::uvec ind;
    int k;
    if (!processed_opts.ind.has_value()) {
        throw std::runtime_error("In glLogisticR, the fields .ind should be specified");
    } else {
        ind = processed_opts.ind.value();

        k = ind.n_elem - 1; // the number of groups
        if (ind(k) != n) {

            throw std::runtime_error("Check opts.ind");
        }
    }

    // Initialize q
    double q;
    if (processed_opts.q == 0) {
        q = 2;
        processed_opts.q = 2;
    } else {
        q = processed_opts.q;
        if (q < 1) {
            throw std::runtime_error("q should be larger than 1");
        }
    }

    // The parameter 'weight' contains the weight for each training sample.
    // The summation of the weights for all the samples equals to 1.
    arma::vec weight(m);
    if (processed_opts.sWeight.has_value()) {
        arma::vec sWeight = processed_opts.sWeight.value();

        if (sWeight.n_elem != 2 || sWeight(0) <= 0 || sWeight(1) <= 0) {
            throw std::runtime_error("Check opts.sWeight, which contains two positive values");
        }

        // we process the weight, so that the summation of the weights for all
        // the samples is 1.

        arma::uvec p_flag = arma::find(y == 1);                  // the indices of the postive samples
        double m1 = p_flag.n_elem * sWeight(0);    // the total weight for the positive samples
        double m2 = (m - p_flag.n_elem) * sWeight(1);   // the total weight for the negative samples

        weight.fill(0);
        weight.elem(p_flag).fill(sWeight(0) / (m1 + m2));
        arma::uvec n_flag = arma::find(y != 1);
        weight.elem(n_flag).fill(sWeight(1) / (m1 + m2));
    } else {
        weight.ones();
        weight /= m;             // if not specified, we apply equal weight
    }
    // gWeight: the weigtht for each group
    arma::vec gWeight(k);
    if (processed_opts.gWeight.has_value()) {
        gWeight = processed_opts.gWeight.value();
        if (gWeight.n_elem != k) {
            throw std::runtime_error("opts.gWeight should a " + std::to_string(k) + " x 1 vector");
        }

        if (arma::min(gWeight) <= 0) {
            throw std::runtime_error(".gWeight should be positive");
        }
    } else {
        gWeight.ones();
    }

    // Starting point initialization

    arma::uvec p_flag = arma::find(y == 1);                  // the indices of the postive samples
    double m1 = arma::sum(weight.elem(p_flag));         // the total weight for the positive samples
    double m2 = 1 - m1;                        // the total weight for the negative samples
    // process the regularization parameter
    double lambda;
    if (processed_opts.rFlag == 0) {
        lambda = z;
    } else { // z here is the scaling factor lying in [0,1]
        //     if (z<0 || z>1)
        //         error('\n opts.rFlag=1, and z should be in [0,1]');
        //     end

        // we compute ATb for computing lambda_max, when the input z is a ratio

        arma::vec b(m);
        b.elem(p_flag).fill(m2);
        arma::uvec n_flag = arma::find(y != 1);
        b.elem(n_flag).fill(-m1);
        b = b % weight;

        // compute AT b
        arma::vec ATb;
        if (processed_opts.nFlag == 0) {
            ATb = A.t() * b;
        } else if (processed_opts.nFlag == 1) {
            ATb = A.t() * b - arma::sum(b) * mu;
            ATb = ATb / nu;
        } else {
            arma::vec invNu = b / nu;
            ATb = A.t() * invNu - arma::sum(invNu) * mu;
        }

        // the conjugate of q
        double q_bar;
        if (q == 1) {
            q_bar = arma::datum::inf;
        } else if (q >= 1e6) {     // when q>=1e6, we treat it as inf
            q_bar = 1;
        } else {
            q_bar = q / (q - 1);
        }

        // compute the norm of ATb corresponding to each group
        arma::vec norm_ATb(k);
        for (int i = 0; i < k; i++) {
            arma::vec group_ATb = ATb.subvec(ind(i), ind(i+1) - 1);
            norm_ATb(i) = arma::norm(group_ATb, q_bar);
        }

        // incorporate the gWeight
        norm_ATb = norm_ATb / gWeight;

        // compute lambda_max
        double lambda_max = arma::max(norm_ATb);

        // As .rFlag=1, we set lambda as a ratio of lambda_max
        lambda = z * lambda_max;
    }
    // initialize a starting point
    arma::vec x(n);
    double c;
    if (processed_opts.init == 2) {
        x.zeros();
        c = std::log(m1 / m2);
    } else {
        if (processed_opts.x0.has_value()) {
            x = processed_opts.x0.value();
            if (x.n_elem != n) {
                throw std::runtime_error("Check the input .x0");
            }
        } else {
            x.zeros();
        }

        if (processed_opts.c0.has_value() && processed_opts.c0 != 0) {
            c = processed_opts.c0.value();
        } else {
            c = std::log(m1 / m2);
        }
    }

    // compute A x
    arma::vec Ax(m);
    if (processed_opts.nFlag == 0) {
        Ax = A * x;
    } else if (processed_opts.nFlag == 1) {
        arma::vec invNu = x / nu;
        double mu_invNu = arma::dot(mu, invNu);
        Ax = A * invNu - mu_invNu;
    } else {
        Ax = A * x - arma::dot(mu, x);
        Ax = Ax / nu;
    }
    // Initialize output vectors
    arma::vec funVal(processed_opts.maxIter);
    arma::vec ValueL(processed_opts.maxIter);
    int iterStep = 0;

    // The main program

    if (processed_opts.mFlag == 0 && processed_opts.lFlag == 0) {

        bool bFlag = false; // this flag tests whether the gradient step only changes a little

        double L = 1.0 / m; // the intial guess of the Lipschitz continuous gradient

        arma::vec weighty = weight % y;
        // the product between weight and y

        // assign xp with x, and Axp with Ax
        arma::vec xp = x, Axp = Ax, xxp(n, arma::fill::zeros);
        double cp = c, ccp = 0;

        // The Armijo Goldstein line search schemes + accelearted gradient descent

        double alphap = 0, alpha = 1;

        for (iterStep = 1; iterStep <= processed_opts.maxIter; iterStep++) {
            // --------------------------- step 1 ---------------------------
            // compute search point s based on xp and x (with beta)
            double beta = (alphap - 1) / alpha;
            arma::vec s = x + beta * xxp;
            double sc = c + beta * ccp;

            // --------------------------- step 2 ---------------------------
            // line search for L and compute the new approximate solution x

            // compute As=A*s
            arma::vec As = Ax + beta * (Ax - Axp);

            // aa= - diag(y) * (A * s + sc)
            arma::vec aa = -y % (As + sc);

            // fun_s is the logistic loss at the search point
            arma::vec bb = arma::max(aa, arma::zeros(aa.n_elem));
            double fun_s = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

            // compute prob=[p_1;p_2;...;p_m]
            arma::vec prob = 1.0 / (1 + arma::exp(aa));

            // b= - diag(y.* weight) * (1 - prob)
            arma::vec b = -weighty % (1 - prob);

            double gc = arma::sum(b); // the gradient of c

            // compute g= AT b, the gradient of x
            arma::vec g;
            if (processed_opts.nFlag == 0) {
                g = A.t() * b;
            } else if (processed_opts.nFlag == 1) {
                g = A.t() * b - arma::sum(b) * mu;
                g = g / nu;
            } else {
                arma::vec invNu = b / nu;
                g = A.t() * invNu - arma::sum(invNu) * mu;
            }

            // copy x and Ax to xp and Axp
            xp = x;    Axp = Ax;
            cp = c;
            double fun_x;
            while (true) {
                // let s walk in a step in the antigradient of s to get v
                // and then do the L1/Lq-norm regularized projection
                arma::vec v = s - g / L;
                c = sc - gc / L;

                // L1/Lq-norm regularized projection
                if (q < 1e6) {
                    x = eppVector(v, ind, k, n, lambda / L * gWeight, q);
                } else { // when q>=1e6, we treat q as inf
                    x = eppVector(v, ind, k, n, lambda / L * gWeight, 1e6);
                }

                v = x - s;  // the difference between the new approximate solution x
                // and the search point s

                // compute A x
                if (processed_opts.nFlag == 0) {
                    Ax = A * x;
                } else if (processed_opts.nFlag == 1) {
                    arma::vec invNu = x / nu;
                    double mu_invNu = arma::dot(mu, invNu);
                    Ax = A * invNu - mu_invNu;
                } else {
                    Ax = A * x - arma::dot(mu, x);
                    Ax = Ax / nu;
                }

                // aa= - diag(y) * (A * x + c)
                aa = -y % (Ax + c);

                // fun_x is the logistic loss at the new approximate solution
                bb = arma::max(aa, arma::zeros(aa.n_elem));
                fun_x = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

                double r_sum = (arma::dot(v, v) + std::pow(c - sc, 2)) / 2;
                double l_sum = fun_x - fun_s - arma::dot(v, g) - (c - sc) * gc;

                if (r_sum <= 1e-20) {
                    bFlag = true; // this shows that, the gradient step makes little improvement
                    break;
                }

                // the condition is fun_x <= fun_s + v'* g + c * gc
                //                           + L/2 * (v'*v + (c-sc)^2 )
                if (l_sum <= r_sum * L) {
                    break;
                } else {
                    L = std::max(2 * L, l_sum / r_sum);
                    //fprintf('\n L=%e, r_sum=%e',L, r_sum);
                }
            }

            // --------------------------- step 3 ---------------------------
            // update alpha and alphap, and check whether converge
            alphap = alpha;
            alpha = (1 + std::sqrt(4 * alpha * alpha + 1)) / 2;

            ValueL(iterStep - 1) = L;
            // store values for L

            xxp = x - xp;    ccp = c - cp;

            // the q-norm of x
            arma::vec norm_x_k(k);
            for (int i = 0; i < k; i++) {
                arma::vec x_group = x.subvec(ind(i), ind(i+1) - 1);
                norm_x_k(i) = arma::norm(x_group, q);
            }

            // function value = loss + regularizatioin
            funVal(iterStep - 1) = fun_x + lambda * arma::dot(norm_x_k, gWeight);

            if (bFlag) {
                // fprintf('\n The program terminates as the gradient step changes the solution very small.');
                break;
            }
            switch (processed_opts.tFlag) {
                case 0:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <= processed_opts.tol) {
                            goto exit_loop_1;
                        }
                    }
                    break;
                case 1:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <=
                            processed_opts.tol * funVal(iterStep - 2)) {
                            goto exit_loop_1;
                        }
                    }
                    break;
                case 2:
                    if (funVal(iterStep - 1) <= processed_opts.tol) {
                        goto exit_loop_1;
                    }
                    break;
                case 3:
                    {
                        double norm_xxp = arma::norm(xxp);
                        if (norm_xxp <= processed_opts.tol) {
                            goto exit_loop_1;
                        }
                    }
                    break;
                case 4:
                    {
                        double norm_xp = arma::norm(xp);
                        double norm_xxp = arma::norm(xxp);
                        if (norm_xxp <= processed_opts.tol * std::max(norm_xp, 1.0)) {
                            goto exit_loop_1;
                        }
                    }
                    break;
                case 5:
                    if (iterStep >= processed_opts.maxIter) {
                        goto exit_loop_1;
                    }
                    break;
            }
        }
        exit_loop_1:;
    }


    // Reformulated problem + Nemirovski's line search scheme


    if (processed_opts.mFlag == 1 && processed_opts.lFlag == 0 && processed_opts.q == 2) {

        bool bFlag = false; // this flag tests whether the gradient step only changes a little

        double L = 1.0 / m; // the intial guess of the Lipschitz continuous gradient

        arma::vec weighty = weight % y;
        // the product between weight and y

        // assign xp with x, and Axp with Ax
        arma::vec xp = x, Axp = Ax, xxp(n, arma::fill::zeros);
        double cp = c, ccp = 0;

        arma::vec t(k);
        for (int i = 0; i < k; i++) {
            arma::vec x_group = x.subvec(ind(i), ind(i+1) - 1);
            t(i) = arma::norm(x_group, 2);
        }
        arma::vec tp = t;
        // t is the upper bound of the 2-norm of x

        double alphap = 0, alpha = 1;

        for (iterStep = 1; iterStep <= processed_opts.maxIter; iterStep++) {
            // --------------------------- step 1 ---------------------------
            // compute search point s based on xp and x (with beta)
            double beta = (alphap - 1) / alpha;
            arma::vec s = x + beta * xxp;
            double sc = c + beta * ccp;
            arma::vec s_t = t + beta * (t - tp);

            // --------------------------- step 2 ---------------------------
            // line search for L and compute the new approximate solution x

            // compute As=A*s
            arma::vec As = Ax + beta * (Ax - Axp);

            // aa= - diag(y) * (A * s + sc)
            arma::vec aa = -y % (As + sc);

            // fun_s is the logistic loss at the search point
            arma::vec bb = arma::max(aa, arma::zeros(aa.n_elem));
            double fun_s = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

            // compute prob=[p_1;p_2;...;p_m]
            arma::vec prob = 1.0 / (1 + arma::exp(aa));

            // b= - diag(y.* weight) * (1 - prob)
            arma::vec b = -weighty % (1 - prob);

            double gc = arma::sum(b); // the gradient of c

            // compute g= AT b, the gradient of x
            arma::vec g;
            if (processed_opts.nFlag == 0) {
                g = A.t() * b;
            } else if (processed_opts.nFlag == 1) {
                g = A.t() * b - arma::sum(b) * mu;
                g = g / nu;
            } else {
                arma::vec invNu = b / nu;
                g = A.t() * invNu - arma::sum(invNu) * mu;
            }

            // copy x and Ax to xp and Axp
            xp = x;    Axp = Ax;
            cp = c;    tp = t;
            double fun_x;

            while (true) {
                // let s walk in a step in the antigradient of s to get v
                // and then do the L1/Lq-norm regularized projection
                c = sc - gc / L;
                arma::vec u = s - g / L;
                arma::vec v = s_t - lambda / L;

                // projection
                auto result = eppVectorR(u, v, ind, n, k);
                x = result.first;
                t = result.second;

                arma::vec v_diff = x - s;  // the difference between the new approximate solution x
                // and the search point s
                arma::vec v_t = t - s_t;

                // compute A x
                if (processed_opts.nFlag == 0) {
                    Ax = A * x;
                } else if (processed_opts.nFlag == 1) {
                    arma::vec invNu = x / nu;
                    double mu_invNu = arma::dot(mu, invNu);
                    Ax = A * invNu - mu_invNu;
                } else {
                    Ax = A * x - arma::dot(mu, x);
                    Ax = Ax / nu;
                }

                // aa= - diag(y) * (A * x + c)
                aa = -y % (Ax + c);

                // fun_x is the logistic loss at the new approximate solution
                bb = arma::max(aa, arma::zeros(aa.n_elem));
                fun_x = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

                double r_sum = (arma::dot(v_diff, v_diff) + std::pow(c - sc, 2) + arma::dot(v_t, v_t)) / 2;
                double l_sum = fun_x - fun_s - arma::dot(v_diff, g) - (c - sc) * gc;

                if (r_sum <= 1e-20) {
                    bFlag = true; // this shows that, the gradient step makes little improvement
                    break;
                }

                // the condition is fun_x <= fun_s + v'* g + c * gc
                //                           + L/2 * (v'*v + (c-sc)^2 )
                if (l_sum <= r_sum * L) {
                    break;
                } else {
                    L = std::max(2 * L, l_sum / r_sum);
                    //fprintf('\n L=%e, r_sum=%e',L, r_sum);
                }
            }

            // --------------------------- step 3 ---------------------------
            // update alpha and alphap, and check whether converge
            alphap = alpha;
            alpha = (1 + std::sqrt(4 * alpha * alpha + 1)) / 2;

            ValueL(iterStep - 1) = L;
            // store values for L

            xxp = x - xp;    ccp = c - cp;
            funVal(iterStep - 1) = fun_x + lambda * arma::dot(t, gWeight);

            if (bFlag) {
                // fprintf('\n The program terminates as the gradient step changes the solution very small.');
                break;
            }

            switch (processed_opts.tFlag) {
                case 0:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <= processed_opts.tol) {
                            goto exit_loop_2;
                        }
                    }
                    break;
                case 1:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <=
                            processed_opts.tol * funVal(iterStep - 2)) {
                            goto exit_loop_2;
                        }
                    }
                    break;
                case 2:
                    if (funVal(iterStep - 1) <= processed_opts.tol) {
                        goto exit_loop_2;
                    }
                    break;
                case 3:
                    {
                        double norm_xxp = std::sqrt(arma::dot(xxp, xxp) + std::pow(arma::norm(t - tp), 2) + std::pow(c - cp, 2));
                        if (norm_xxp <= processed_opts.tol) {
                            goto exit_loop_2;
                        }
                    }
                    break;
                case 4:
                    {
                        double norm_xp = std::sqrt(arma::dot(xp, xp) + cp * cp + arma::dot(tp, tp));
                        double norm_xxp = std::sqrt(arma::dot(xxp, xxp) + std::pow(arma::norm(t - tp), 2) + std::pow(c - cp, 2));
                        if (norm_xxp <= processed_opts.tol * std::max(norm_xp, 1.0)) {
                            goto exit_loop_2;
                        }
                    }
                    break;
                case 5:
                    if (iterStep >= processed_opts.maxIter) {
                        goto exit_loop_2;
                    }
                    break;
            }
        }
        exit_loop_2:;
    }


    // Reformulated problem + adaptive line search scheme

    if (processed_opts.mFlag == 1 && processed_opts.lFlag == 1 && processed_opts.q == 2) {

        bool bFlag = false; // this flag tests whether the gradient step only changes a little

        double L = 1.0 / m; // the intial guess of the Lipschitz continuous gradient

        arma::vec weighty = weight % y;
        // the product between weight and y

        double gamma = 1;
        // we shall set the value of gamma = L,
        // and L is appropriate for the starting point

        // assign xp with x, and Axp with Ax
        arma::vec xp = x, Axp = Ax, xxp(n, arma::fill::zeros);
        double cp = c, ccp = 0;

        arma::vec t(k);
        for (int i = 0; i < k; i++) {
            arma::vec x_group = x.subvec(ind(i), ind(i+1) - 1);
            t(i) = arma::norm(x_group, 2);
        }
        arma::vec tp = t;
        // t is the upper bound of the 2-norm of x

        double alphap = 0;

        for (iterStep = 1; iterStep <= processed_opts.maxIter; iterStep++) {
            double alpha, beta;
            arma::vec s, s_t, As;
            double sc;
            double fun_x;

            while (true) {
                if (iterStep != 1) {
                    alpha = (-gamma + std::sqrt(gamma * gamma + 4 * L * gamma)) / (2 * L);
                    beta = (gamma - gamma * alphap) / (alphap * gamma + alphap * L * alpha);
                    // beta is the coefficient for generating search point s

                    s = x + beta * xxp;
                    s_t = t + beta * (t - tp);
                    sc = c + beta * ccp;
                    As = Ax + beta * (Ax - Axp);
                } else {
                    alpha = (-1 + std::sqrt(5)) / 2;
                    beta = 0;
                    s = x;
                    s_t = t;
                    sc = c;
                    As = Ax;
                }

                // aa= - diag(y) * (A * s + sc)
                arma::vec aa = -y % (As + sc);

                // fun_s is the logistic loss at the search point
                arma::vec bb = arma::max(aa, arma::zeros(aa.n_elem));
                double fun_s = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

                // compute prob=[p_1;p_2;...;p_m]
                arma::vec prob = 1.0 / (1 + arma::exp(aa));

                // b= - diag(y.* weight) * (1 - prob)
                arma::vec b = -weighty % (1 - prob);

                double gc = arma::sum(b); // the gradient of c

                // compute g= AT b, the gradient of x
                arma::vec g;
                if (processed_opts.nFlag == 0) {
                    g = A.t() * b;
                } else if (processed_opts.nFlag == 1) {
                    g = A.t() * b - arma::sum(b) * mu;
                    g = g / nu;
                } else {
                    arma::vec invNu = b / nu;
                    g = A.t() * invNu - arma::sum(invNu) * mu;
                }

                // let s walk in a step in the antigradient of s to get v
                // and then do the L1/Lq-norm regularized projection
                double cnew = sc - gc / L;
                arma::vec u = s - g / L;
                arma::vec v = s_t - lambda / L;

                // projection
                auto result = eppVectorR(u, v, ind, n, k);
                arma::vec xnew = result.first;
                arma::vec tnew = result.second;

                arma::vec v_diff = xnew - s;  // the difference between the new approximate solution x
                // and the search point s
                arma::vec v_t = tnew - s_t;

                // compute A xnew
                arma::vec Axnew;
                if (processed_opts.nFlag == 0) {
                    Axnew = A * xnew;
                } else if (processed_opts.nFlag == 1) {
                    arma::vec invNu = xnew / nu;
                    double mu_invNu = arma::dot(mu, invNu);
                    Axnew = A * invNu - mu_invNu;
                } else {
                    Axnew = A * xnew - arma::dot(mu, xnew);
                    Axnew = Axnew / nu;
                }

                // aa= - diag(y) * (A * x + c)
                aa = -y % (Axnew + cnew);

                // fun_x is the logistic loss at the new approximate solution
                bb = arma::max(aa, arma::zeros(aa.n_elem));
                fun_x = arma::dot(weight, arma::log(arma::exp(-bb) + arma::exp(aa - bb)) + bb);

                double r_sum = (arma::dot(v_diff, v_diff) + std::pow(cnew - sc, 2) + arma::dot(v_t, v_t)) / 2;
                double l_sum = fun_x - fun_s - arma::dot(v_diff, g) - (cnew - sc) * gc;

                if (r_sum <= 1e-20) {
                    bFlag = true; // this shows that, the gradient step makes little improvement
                    break;
                }

                // the condition is fun_x <= fun_s + v'* g + c * gc
                //                           + L/2 * (v'*v + (c-sc)^2 )
                if (l_sum <= r_sum * L) {
                    // Update all variables here when condition is satisfied
                    gamma = L * alpha * alpha;
                    alphap = alpha;
                    // update gamma, and alphap

                    double tao = L * r_sum / l_sum;
                    if (tao >= 5) {
                        L = L * 0.8;
                    }
                    // decrease the value of L

                    xp = x;  x = xnew; xxp = x - xp;
                    Axp = Ax; Ax = Axnew;
                    tp = t; t = tnew;
                    cp = c; c = cnew; ccp = c - cp;
                    // update x, t, c, and Ax

                    break;
                } else {
                    L = std::max(2 * L, l_sum / r_sum);
                    //fprintf('\n L=%e, r_sum=%e',L, r_sum);
                }
            }

            ValueL(iterStep - 1) = L;
            // store values for L

            funVal(iterStep - 1) = fun_x + lambda * arma::dot(t, gWeight);

            if (bFlag) {
                // fprintf('\n The program terminates as the gradient step changes the solution very small.');
                break;
            }

            switch (processed_opts.tFlag) {
                case 0:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <= processed_opts.tol) {
                            goto exit_loop_3;
                        }
                    }
                    break;
                case 1:
                    if (iterStep >= 2) {
                        if (std::abs(funVal(iterStep - 1) - funVal(iterStep - 2)) <=
                            processed_opts.tol * funVal(iterStep - 2)) {
                            goto exit_loop_3;
                        }
                    }
                    break;
                case 2:
                    if (funVal(iterStep - 1) <= processed_opts.tol) {
                        goto exit_loop_3;
                    }
                    break;
                case 3:
                    {
                        double norm_xxp = std::sqrt(arma::dot(xxp, xxp) + std::pow(arma::norm(t - tp), 2) + std::pow(c - cp, 2));
                        if (norm_xxp <= processed_opts.tol) {
                            goto exit_loop_3;
                        }
                    }
                    break;
                case 4:
                    {
                        double norm_xp = std::sqrt(arma::dot(xp, xp) + cp * cp + arma::dot(tp, tp));
                        double norm_xxp = std::sqrt(arma::dot(xxp, xxp) + std::pow(arma::norm(t - tp), 2) + std::pow(c - cp, 2));
                        if (norm_xxp <= processed_opts.tol * std::max(norm_xp, 1.0)) {
                            goto exit_loop_3;
                        }
                    }
                    break;
                case 5:
                    if (iterStep >= processed_opts.maxIter) {
                        goto exit_loop_3;
                    }
                    break;
            }
        }
        exit_loop_3:;
    }

    if (processed_opts.mFlag == 0 && processed_opts.lFlag == 1) {
        throw std::runtime_error("The function does not support opts.mFlag=0 & opts.lFlag=1!");
    }

    // Resize output vectors to actual number of iterations
    funVal.resize(iterStep);
    ValueL.resize(iterStep);

    glLogisticR_result result;
    result.x = x;
    result.c = c;
    result.funVal = funVal;
    result.ValueL = ValueL;

    return result;
}

