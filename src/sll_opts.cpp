
#include <armadillo>
#include <cmath>
#include <string>
#include <optional>

// Structure to hold all options for Sparse Learning Library
struct sll_opts {
    // Starting point
    std::optional<arma::vec> x0;           // Starting point of x
    std::optional<double> c0;              // Starting point for intercept c
    int init = 2;                          // Initialization method (default: 2)

    // Termination
    int maxIter = 100;                     // Maximum number of iterations (default: 1e4)
    double tol = 1e-3;                     // Tolerance parameter (default: 1e-4, but code uses 1e-3)
    int tFlag = 5;                         // Flag for termination (default: 0)

    // Normalization
    int nFlag = 0;                         // Flag for implicit normalization of A (default: 0)
    std::optional<arma::colvec> mu;        // Row vector to be subtracted from each sample
    std::optional<arma::vec> nu;           // Weight vector for normalization

    // Regularization
    int rFlag = 1;                         // Flag for regularization (default: 0)
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

// Function to validate and set default options for Sparse Learning Library
sll_opts default_sll_opts(sll_opts opts) {
    // Options for Sparse Learning Library
    //
    // Notice:
    // If one or several (even all) fields are empty, sll_opts shall assign the
    // default settings.
    //
    // If some fields of opts have been defined, sll_opts shall check the fields
    // for possible errors.
    //
    //
    // Table of Options.  * * indicates default value.
    //
    // FIELD            DESCRIPTION
    // Starting point
    //
    // .x0               Starting point of x.
    //                   Initialized according to .init.
    //
    // .c0               Starting point for the intercept c (for Logistic Loss)
    //                   Initialized according to .init.
    //
    // .init             .init specifies how to initialize x.
    //                       * 0 => .x0 is set by the function initFactor *
    //                         1 => .x0 and .c0 are defined
    //                         2 => .x0= zeros(n,1), .c0=0
    //
    // Termination
    //
    // .maxIter          Maximum number of iterations.
    //                       *1e4*
    //
    // .tol              Tolerance parameter.
    //                       *1e-4*
    //
    // .tFlag            Flag for termination.
    //                       * 0 => abs( funVal(i)- funVal(i-1) ) <= .tol *
    //                         1 => abs( funVal(i)- funVal(i-1) )
    //                              <= .tol max( funVal(i-1), 1)
    //                         2 => funVal(i) <= .tol
    //                         3 => norm( x_i - x_{i-1}, 2) <= .tol
    //                         4 => norm( x_i - x_{i-1}, 2) <=
    //                              <= .tol max( norm( x_{i-1}, 2), 1 )
    //                         5 => Run the code for .maxIter iterations
    //
    // Normalization
    //
    // .nFlag            Flag for implicit normalization of A.
    //                       * 0 => Do not normalize A *
    //                         1 => A=(A-repmat(mu, m, 1))*diag(nu)^{-1}
    //                         2 => A=diag(nu)^{-1}*(A-repmat(mu,m,1)
    //
    // .mu               Row vector to be substracted from each sample.
    //                           (.mu is used when .nFlag=1 or 2)
    //                       If .mu is not specified, then
    //                            * .mu=mean(A,1) *
    //
    // .nu               Weight (column) vector for normalization
    //                           (.mu is used when .nFlag=1 or 2)
    //                       If .nu is not specified, then
    //                       * .nFlag=1 => .nu=(sum(A.^2, 1)'/m.^{0.5} *
    //                       * .nFlag=2 => .nu=(sum(A.^2, 2)/n.^{0.5} *
    //
    // Regularization
    //
    // .rFlag            Flag for regularization
    //                           (.rFlag is used for the functions with "R")
    //                        * 0 => lambda is the regularization parameter *
    //                          1 => lambda = lambda * lambda_{max}
    //                               where lambda_{max} is the maximum lambda
    //                               that yields the zero solution
    // .rsL2              Regularization parameter value of the squared L2 norm
    //                           (.rsL2 is used only for l1 regularization)
    //                        *.rsL2=0*
    //                    If .rFlag=0, .rsL2 is used without scaling
    //                       .rFlag=1, .rsL2=.rsL2 * lambda_{max}
    //
    // Method & Line Search
    // .lFlag
    //
    // Grooup & Others
    //
    // .ind              Indices for k groups (a k+1 row vector)
    //                   For group lasso only
    //                   Indices for the i-th group are (ind(i)+1):ind(i+1)
    //
    // .q                Value of q in L1/Lq regularization
    //                      *.q=2*
    //
    // .sWeight          The sample (positive and negative) weight
    //                   For the Logistic Loss only
    //                   Positive sample: .sWeight(1)
    //                   Negative sample: sWeight(2)
    //                   *1/m for both positive and negative samples*
    //
    // .gWeight          The weight for different groups
    //                      *.gWeight=1*
    //
    // .fName            The name of the function
    //
    // Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
    //
    // You are suggested to first read the Manual.
    //
    // For any problem, please contact with Jun Liu via j.liu@asu.edu
    //
    // Last modified 7 August 2009.

    // Starting point
    if (opts.init != 0 && opts.init != 1 && opts.init != 2) {
        opts.init = 0; // if .init is not 0, 1, or 2, then use the default 0
    }

    if (!opts.x0.has_value() && opts.init == 1) {
        opts.init = 0; // if .x0 is not defined and .init=1, set .init=0
    }

    // Termination
    if (opts.maxIter < 1) {
        opts.maxIter = 10000;
    }

    // Note: tol is already initialized to 1e-3 in the struct definition

    if (opts.tFlag < 0) {
        opts.tFlag = 0;
    } else if (opts.tFlag > 5) {
        opts.tFlag = 5;
    } else {
        opts.tFlag = static_cast<int>(std::floor(opts.tFlag));
    }

    // Normalization
    if (opts.nFlag != 1 && opts.nFlag != 2) {
        opts.nFlag = 0;
    }

    // Regularization
    if (opts.rFlag != 1) {
        opts.rFlag = 0;
    }

    // Method (Line Search)
    if (opts.lFlag != 1) {
        opts.lFlag = 0;
    }

    if (opts.mFlag != 1) {
        opts.mFlag = 0;
    }

    return opts;
}

