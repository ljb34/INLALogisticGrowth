#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector mu_sparse_cpp(
    double growth,
    double carry_cap,
    double move_const,
    double timestep,
    NumericVector linpoint,
    NumericMatrix grad,
    NumericVector prior_mean,
    int ns,
    int nt,
    const SparseMatrix<double>& CinvG
)
{
    const int N = ns * nt;

    //-----------------------------
    // Build vector a
    //-----------------------------
    VectorXd a(N);
    for (int i = 0; i < N; i++)
        a(i) = growth * std::exp(linpoint[i]) / carry_cap;

    for (int i = 0; i < ns; i++)
        a(i) = 1.0;

    //-----------------------------
    // Build r vector
    //-----------------------------
    VectorXd r(N);

    // first ns entries: prior mean
    for (int i = 0; i < ns; i++)
        r(i) = prior_mean[i];

    // remaining entries: r.vector
    for (int i = ns; i < N; i++) {
        double grad_sq = 0.0;
        for (int j = 0; j < grad.ncol(); j++)
            grad_sq += grad(i, j) * grad(i, j);

        r(i) =
            growth * std::exp(linpoint[i]) * (linpoint[i] - 1.0) / carry_cap
            + growth
            - move_const * grad_sq;
    }

    //-----------------------------
    // Build main spatial block
    //-----------------------------
    SparseMatrix<double> I(ns, ns);
    I.setIdentity();

    SparseMatrix<double> MainBlock =
        (1.0 / timestep) * I + move_const * CinvG;

    //-----------------------------
    // Assemble sparse L
    //-----------------------------
    std::vector< Triplet<double> > triplets;
    triplets.reserve(10 * N);

    // Diagonal (a)
    for (int i = 0; i < N; i++)
        triplets.emplace_back(i, i, a(i));

    // Main blocks (t >= 1)
    for (int t = 1; t < nt; t++) {
        int offset = t * ns;

        for (int k = 0; k < MainBlock.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(MainBlock, k); it; ++it)
                triplets.emplace_back(
                    offset + it.row(),
                    offset + it.col(),
                    it.value()
                );
    }

    // Subdiagonal blocks
    for (int t = 1; t < nt; t++) {
        for (int i = 0; i < ns; i++) {
            triplets.emplace_back(
                t * ns + i,
                (t - 1) * ns + i,
                -1.0 / timestep
            );
        }
    }

    SparseMatrix<double> L(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());
    L.makeCompressed();

    //-----------------------------
    // Solve L * mu = r
    //-----------------------------
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    solver.compute(L);

    if (solver.info() != Success) {
        stop("SparseQR decomposition failed");
    }

    VectorXd mu = solver.solve(r);

    if (solver.info() != Success) {
        stop("Sparse solve failed");
    }

    return wrap(mu);
}