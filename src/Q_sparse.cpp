#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> Q_sparse_cpp(
    double growth,
    double carry_cap,
    double move_const,
    double sigma,
    double timestep,
    NumericVector linpoint,
    int ns,
    int nt,
    const Eigen::SparseMatrix<double>& C,
    const Eigen::SparseMatrix<double>& G,
    const Eigen::SparseMatrix<double>& CinvG,
    const Eigen::SparseMatrix<double>& prior_precision
)
{
    const int N = ns * nt;

    //----------------------------------
    // Precompute Qblock = sigma^2 / h * (C + g G)
    //----------------------------------
    Eigen::SparseMatrix<double> Qblock = C + move_const * G;
    Qblock *= (sigma * sigma / timestep);

    //----------------------------------
    // Precompute main spatial block:
    // (1/h) I + g CinvG
    //----------------------------------
    Eigen::SparseMatrix<double> I(ns, ns);
    I.setIdentity();

    Eigen::SparseMatrix<double> MainBlock =
        (1.0 / timestep) * I + move_const * CinvG;

    //----------------------------------
    // Precompute vector a
    //----------------------------------
    VectorXd a(N);
    for (int i = 0; i < N; i++) {
        a(i) = growth * std::exp(linpoint[i]) / carry_cap;
    }
    for (int i = 0; i < ns; i++) {
        a(i) = 1.0;
    }

    //----------------------------------
    // Triplet accumulator for Q
    //----------------------------------
    std::vector< Triplet<double> > triplets;
    triplets.reserve(10 * N);  // heuristic

    //----------------------------------
    // Loop over time blocks
    //----------------------------------
    for (int t = 0; t < nt; t++) {

        int row_offset = t * ns;

        //----------------------------------
        // Diagonal block: L_tᵀ Qblock L_t
        //----------------------------------
        if (t > 0) {

            // L_t = MainBlock + diag(a_t)
            Eigen::SparseMatrix<double> Lt = MainBlock;

            for (int i = 0; i < ns; i++) {
                Lt.coeffRef(i, i) += a(row_offset + i);
            }

            Eigen::SparseMatrix<double> temp = Qblock * Lt;
            Eigen::SparseMatrix<double> Qt = Lt.transpose() * temp;

            for (int k = 0; k < Qt.outerSize(); ++k)
                for (Eigen::SparseMatrix<double>::InnerIterator it(Qt, k); it; ++it)
                    triplets.emplace_back(
                        row_offset + it.row(),
                        row_offset + it.col(),
                        it.value()
                    );
        }
        else {
            // t = 0: prior precision only
            for (int k = 0; k < prior_precision.outerSize(); ++k)
                for (Eigen::SparseMatrix<double>::InnerIterator it(prior_precision, k); it; ++it)
                    triplets.emplace_back(
                        it.row(),
                        it.col(),
                        it.value()
                    );
        }

        //----------------------------------
        // Off-diagonal coupling with t-1
        //----------------------------------
        if (t > 0) {

            // Subdiagonal: (-1/h) I
            Eigen::SparseMatrix<double> Sub(ns, ns);
            Sub.setIdentity();
            Sub *= (-1.0 / timestep);

            Eigen::SparseMatrix<double> temp = Qblock * Sub;

            for (int k = 0; k < temp.outerSize(); ++k)
                for (Eigen::SparseMatrix<double>::InnerIterator it(temp, k); it; ++it) {

                    int i = row_offset + it.row();
                    int j = row_offset - ns + it.col();

                    triplets.emplace_back(i, j, it.value());
                    triplets.emplace_back(j, i, it.value()); // symmetry
                }
        }
    }

    //----------------------------------
    // Assemble sparse Q
    //----------------------------------
    Eigen::SparseMatrix<double> Q(N, N);
    Q.setFromTriplets(triplets.begin(), triplets.end());
    Q.makeCompressed();

    return Q;
}