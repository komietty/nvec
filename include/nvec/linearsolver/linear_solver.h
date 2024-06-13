#ifndef NVEC_LINEAR_SOLVER_H
#define NVEC_LINEAR_SOLVER_H
#include <iostream>
#include "typedefs.h"
#include <Eigen/Sparse>
#ifdef GC_HAVE_SUITESPARSE
#include "suitesparse_pdefinite.h"
#include "suitesparse_square.h"
#endif

namespace pddg {

inline VecXc solveSquare(SprsC& lhs, VecXc& rhs) {
#ifdef GC_HAVE_SUITESPARSE
    VecXc x;
    auto solver = std::make_unique<solver::SquareSolver<complex>>(lhs);
    solver->solve(x, rhs);
    return x;
#else
    Eigen::SparseLU<SprsC> solver;
    solver.compute(lhs);
    return solver.solve(rhs);
#endif
}

inline VecXc solveSmallestEig(SprsC& L, SprsC& M, int nIter = 50) {
#ifdef GC_HAVE_SUITESPARSE
    auto solver = std::make_unique<pddg::solver::PositiveDefiniteSolver<complex>>(L);
    VecXc u = VecXc::Random(L.rows());
    VecXc x = u;
    for (size_t i = 0; i < nIter; i++) {
        solver->solve(x, M * u);
        double s = std::sqrt(std::abs(x.dot(M * x)));
        x /= s;
        u = x;
    }
    return x;
#else
    Eigen::SparseLU<SprsC> solver;
    solver.compute(L);
    VecXc u = VecXc::Random(L.rows());
    VecXc x = u;
    for (size_t i = 0; i < nIter; i++) {
        x = solver.solve(M * u);
        double s = std::sqrt(std::abs(x.dot(M * x)));
        x /= s;
        u = x;
    }
    return x;
#endif
}

inline VecXd solveSmallestEig(const SprsD& L, const SprsD& M, int nIter = 50) {
    Eigen::SimplicialLDLT<SprsD> solver;
    solver.compute(L);
    srand((unsigned int) 13);
    VecXd u = VecXd::Random(L.rows());
    VecXd x = u;
    for (size_t i = 0; i < nIter; i++) {
        std::cout << "iter " << i << std::endl;
        x = solver.solve(M * u);
        //double s = std::sqrt(std::abs(x.dot(M * x)));
        //x /= s;
        //u = x;
        u = x / std::sqrt(x.dot(M * x));
    }
    return x;
}
}
#endif
