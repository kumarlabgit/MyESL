
#include <armadillo>
#include <string>
#include <iostream>

double initFactor(double x_norm, const arma::mat& Ax, const arma::mat& y, double z,
                  const std::string& funName, double rsL2, double x_2norm) {
    //
    // function initFactor
    //     compute the an optimal constant factor for the initialization
    //
    //
    // Input parameters:
    // x_norm-      the norm of the starting point
    // Ax-          A*x, with x being the initialization point
    // y-           the response matrix
    // z-           the regularization parameter or the ball
    // funName-     the name of the function
    //
    // Output parameter:
    // ratio-       the computed optimal initialization point is ratio*x
    //
    // Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
    //
    // For any problem, please contact with Jun Liu via j.liu@asu.edu
    //
    // Last revised on August 2, 2009.

    double ratio = 0.0;

    if (funName == "LeastC") {
        double ratio_max = z / x_norm;
        double ratio_optimal = arma::as_scalar(Ax.t() * y) / (arma::as_scalar(Ax.t() * Ax) + rsL2 * x_2norm);

        if (std::abs(ratio_optimal) <= ratio_max) {
            ratio = ratio_optimal;
        } else if (ratio_optimal < 0) {
            ratio = -ratio_max;
        } else {
            ratio = ratio_max;
        }
        // std::cout << "\n ratio=" << ratio << "," << ratio_optimal << "," << ratio_max << std::endl;

    } else if (funName == "LeastR") {
        ratio = (arma::as_scalar(Ax.t() * y) - z * x_norm) / (arma::as_scalar(Ax.t() * Ax) + rsL2 * x_2norm);
        // std::cout << "\n ratio=" << ratio << std::endl;

    } else if (funName == "glLeastR") {
        ratio = (arma::as_scalar(Ax.t() * y) - z * x_norm) / arma::as_scalar(Ax.t() * Ax);
        // std::cout << "\n ratio=" << ratio << std::endl;

    } else if (funName == "mcLeastR") {
        arma::vec Ax_vec = arma::vectorise(Ax);
        arma::vec y_vec = arma::vectorise(y);
        ratio = (arma::as_scalar(Ax_vec.t() * y_vec) - z * x_norm) / std::pow(arma::norm(Ax, "fro"), 2);
        // std::cout << "\n ratio=" << ratio << std::endl;

    } else if (funName == "mtLeastR") {
        ratio = (arma::as_scalar(Ax.t() * y) - z * x_norm) / arma::as_scalar(Ax.t() * Ax);
        // std::cout << "\n ratio=" << ratio << std::endl;

    } else if (funName == "nnLeastR") {
        ratio = (arma::as_scalar(Ax.t() * y) - z * x_norm) / (arma::as_scalar(Ax.t() * Ax) + rsL2 * x_2norm);
        ratio = std::max(0.0, ratio);

    } else if (funName == "nnLeastC") {
        double ratio_max = z / x_norm;
        double ratio_optimal = arma::as_scalar(Ax.t() * y) / (arma::as_scalar(Ax.t() * Ax) + rsL2 * x_2norm);

        if (ratio_optimal < 0) {
            ratio = 0;
        } else if (ratio_optimal <= ratio_max) {
            ratio = ratio_optimal;
        } else {
            ratio = ratio_max;
        }
        // std::cout << "\n ratio=" << ratio << "," << ratio_optimal << "," << ratio_max << std::endl;

    } else if (funName == "mcLeastC") {
        double ratio_max = z / x_norm;
        arma::vec Ax_vec = arma::vectorise(Ax);
        arma::vec y_vec = arma::vectorise(y);
        double ratio_optimal = arma::as_scalar(Ax_vec.t() * y_vec) / std::pow(arma::norm(Ax.t() * Ax, "fro"), 2);

        if (std::abs(ratio_optimal) <= ratio_max) {
            ratio = ratio_optimal;
        } else if (ratio_optimal < 0) {
            ratio = -ratio_max;
        } else {
            ratio = ratio_max;
        }

    } else {
        std::cout << "\n The specified funName is not supprted" << std::endl;
    }

    return ratio;
}

