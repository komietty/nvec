#include "nvec/vert_rosy_field.h"
#include "nvec/linearsolver/linear_solver.h"
#include <igl/readOBJ.h>
#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
using namespace pddg;

void compute_cotan_laplacian(const Hmesh& G, SprsD & L) {
    std::vector<TripD> T;
    for (size_t iE = 0; iE < G.nE; iE++) {
        auto iH  = G.edge2half[iE];
        auto iVt = G.tail[iH];
        auto iVh = G.head[iH];
        auto w = G.edgeCotan[iE];
        T.emplace_back(iVt, iVt,  w);
        T.emplace_back(iVh, iVh,  w);
        T.emplace_back(iVt, iVh, -w);
        T.emplace_back(iVh, iVt, -w);
    }
    L.resize(G.nV, G.nV);
    L.setFromTriplets(T.begin(), T.end());
}

VecXd compute(const Hmesh& mesh, const VecXd& rho){
    SprsD L;
    SprsD M;
    compute_cotan_laplacian(mesh, L);
    M = static_cast<SprsD>(mesh.baryDualArea.asDiagonal());
    double a = mesh.faceArea.sum();
    double r = (M * rho).sum();
    VecXd bar = VecXd::Constant(M.rows(), r / a);
    Eigen::SimplicialLDLT<SprsD> llt;
    llt.compute(L);
    return llt.solve(-M * (rho - bar));
}

inline void computePoisson(const Hmesh& mesh, MatXd& field, VecXd& oneform) {
    oneform.resize(mesh.nE);
    VecXd rho1(mesh.nV);
    VecXd rho2(mesh.nV);
    rho1.setZero();
    rho2.setZero();
    rho1[0] = 1e3;
    rho2[100] = 1e3;
    field.resize(mesh.nF, 3);
    VecXd scalarPot = compute(mesh, rho1);
    VecXd vectorPot = compute(mesh, rho2);
    for(Face f: mesh.faces) {
        auto a = f.area();
        Vec3d n = f.normal();
        Vec3d o {0, 0, 0};

        for (Half h: f.adjHalfs()) {
            Vert v = h.prev().tail();
            Row3d v1 = h.vec();
            Row3d v2 = n.cross(v1);
            o += v2 * scalarPot[v.id] / (a * 2);
            o += v1 * vectorPot[v.id] / (a * 2);
        }
        Row3d c = f.center();
        Row3d u {-c.z(), 0, c.x()};
        u -= n * u.dot(n);
        o += u.normalized() * 0.3;
        field.row(f.id) = o.normalized();
    }

    for (Edge e: mesh.edges) {
        Half h  = e.half();
        Face fa = h.face();
        Face fb = h.twin().face();
        Row3d f1 {0, 0 ,0};
        Row3d f2 {0, 0 ,0};
        if (!h.isBoundary())        f1 = field.row(fa.id);
        if (!h.twin().isBoundary()) f2 = field.row(fb.id);
        oneform[e.id] = (f1 + f2).dot(h.vec());
    }
}

void computeHodgeDecomposition(
        const Hmesh& mesh,
        const VecXd& src,
        VecXd& E,
        VecXd& C,
        VecXd& H
) {
    E.resize(mesh.nE);
    C.resize(mesh.nE);
    H.resize(mesh.nE);
    SprsD d0 = mesh.d0();
    SprsD d1 = mesh.d1();
    SprsD h1 = mesh.h1();
    SprsD h1i = mesh.h1i();
    SprsD lhsA = d0.transpose() * mesh.h1() * d0;
    SprsD lhsB = d1 * h1i * d1.transpose();
    SprsD I(mesh.nV, mesh.nV);
    I.setIdentity();
    Eigen::SimplicialLLT<SprsD> chol;
    Eigen::SparseLU<SprsD> lu;
    chol.compute(lhsA + I * 1e-10);
    lu.compute(lhsB);
    E = d0 * chol.solve(d0.transpose() * h1 * src);
    C = h1i * d1.transpose() * lu.solve(d1 * src);
    H = src - E - C;
}

double pointwiseCurvature(const Hmesh& m, Edge e) {
    Vert v1 = e.half().tail();
    Vert v2 = e.half().head();
    return (m.angleDefect(v1.id) / v1.circArea() +
            m.angleDefect(v2.id) / v2.circArea()) * 0.5;
}

VecXd computeKillingField(const Hmsh& mesh) {
    SprsD h0  = mesh.h0();
    SprsD h1  = mesh.h1(); // B
    SprsD h2  = mesh.h2();
    SprsD d0  = mesh.d0();
    SprsD d1  = mesh.d1();
    SprsD h0i = mesh.h0i();
    SprsD h1i = mesh.h1i();
    SprsD d0t = d0.transpose();
    SprsD d1t = d1.transpose();
    SprsD adjoint1  = h0i * d0t * h1;
    SprsD adjoint2  = h1i * d1t * h2;
    SprsD laplacian = adjoint2 * d1 + d0 * adjoint1;
    VecXd tmp(mesh.nE);
    for (Edge e: mesh.edges)
        tmp[e.id] = pointwiseCurvature(mesh, e);
    SprsD diag = static_cast<SprsD>(VecXd::Constant(h1.rows(), -1e-10).asDiagonal());
    SprsD G = static_cast<SprsD>(tmp.asDiagonal());
    SprsD R = laplacian + d0 * adjoint1 - 2. * h1 * G;
    std::cout << (R.toDense() - R.transpose()).norm() << std::endl;
    /*
    MatXd temp1(6, 6);
    MatXd temp2(6, 6);
    temp1 <<
    4.95576, 1.33333, 1.33333, -1.33333, -1.33334,        0,
    1.33333, 4.95575, 1.33333,        0,  1.33334, -1.33333,
    1.33333, 1.33333, 4.95576,  1.33333,        0,  1.33333,
    -1.33333,      0, 1.33333,  4.95576,  1.33334,  1.33333,
    -1.33334, 1.33334,      0,  1.33334,  4.95576, -1.33334,
    0,      -1.33333, 1.33333,  1.33333, -1.33334,  4.95575;

    temp2 <<
    0.57735,        0,        0,       0,       0,        0,
    0,       0.577351,        0,       0,       0,        0,
    0,              0, 0.577351,       0,       0,        0,
    0,              0,        0, 0.57735,       0,        0,
    0,              0,        0,       0, 0.57735,        0,
    0,              0,        0,       0,       0, 0.577351;

    SprsD tempS1 = temp1.transpose().sparseView();
    SprsD tempS2 = temp2.transpose().sparseView();
    std::cout << solveSmallestEig(tempS1, tempS2, 200) << std::endl;
    std::cout << "..." << std::endl;
    std::cout << solveSmallestEig(tempS1, tempS2, 200) << std::endl;
    std::cout << "..." << std::endl;
    std::cout << solveSmallestEig(tempS1, tempS2, 200) << std::endl;
    */
    return solveSmallestEig(R, h1, 20);
}

MatXd computeWhitenyInterporaiton(const Hmsh& g, const VecXd& oneForm) {
    MatXd F(g.nF, 3);
    for (Face f: g.faces) {
        Half cH = f.half();
        Half nH = cH.next();
        Half pH = cH.prev();
        Row3d pi = cH.tail().pos();
        Row3d pj = cH.head().pos();
        Row3d pk = nH.head().pos();
        Row3d eij = pj - pi;
        Row3d ejk = pk - pj;
        Row3d eki = pi - pk;
        double cij = oneForm[cH.edge().id] * (cH.isCanonical() ? 1. : -1.);
        double cjk = oneForm[nH.edge().id] * (nH.isCanonical() ? 1. : -1.);
        double cki = oneForm[pH.edge().id] * (pH.isCanonical() ? 1. : -1.);
        Row3d a = (eki - ejk) * cij;
        Row3d b = (eij - eki) * cjk;
        Row3d c = (ejk - eij) * cki;
        F.row(f.id) = f.normal().cross(a + b + c) / (6 * f.area());
    }
    return F;
}

int main(int argc, char *argv[]) {
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    MatXd V;
    MatXi F;
    igl::readOBJ("/Users/komietty/dev/models/spot.obj", V, F);
    auto mesh = std::make_unique<Hmsh>(V, F);

    polyscope::init();
    auto surf = polyscope::registerSurfaceMesh("mesh", mesh->pos, mesh->idx);

    VecXd killing = computeKillingField(*mesh);

    //MatXd ambient;
    //VecXd oneform;
    //VecXd exact;
    //VecXd coext;
    //VecXd harmo;
    //computePoisson(*mesh, ambient, oneform);
    //computeHodgeDecomposition(*mesh, oneform, exact, coext, harmo);
    //Row3d v_ex = computeWhitenyInterporaiton(*mesh, exact).colwise().mean();
    //Row3d v_cx = computeWhitenyInterporaiton(*mesh, coext).colwise().mean();
    //Row3d v_hm = computeWhitenyInterporaiton(*mesh, harmo).colwise().mean();
    //std::cout << v_ex.norm() << std::endl;
    //std::cout << v_cx.norm() << std::endl;
    //std::cout << v_hm.norm() << std::endl;
    //auto tq1 = surf->addFaceVectorQuantity("ambient 1", ambient);
    //auto tq2 = surf->addFaceVectorQuantity("ambient 2", computeWhitenyInterporaiton(*mesh, oneform));
    //auto tq3 = surf->addFaceVectorQuantity("exact", computeWhitenyInterporaiton(*mesh, exact));
    //auto tq4 = surf->addFaceVectorQuantity("coext", computeWhitenyInterporaiton(*mesh, coext));
    //auto tq5 = surf->addFaceVectorQuantity("harmo", computeWhitenyInterporaiton(*mesh, harmo));
    auto tq6 = surf->addFaceVectorQuantity("killing", computeWhitenyInterporaiton(*mesh, killing));
    tq6->setEnabled(true);
    //tq1->setVectorLengthScale(0.02);
    //tq2->setVectorLengthScale(0.02);
    //tq3->setVectorLengthScale(0.02);
    //tq4->setVectorLengthScale(0.01);
    //tq5->setVectorLengthScale(0.02);

    polyscope::show();
}
