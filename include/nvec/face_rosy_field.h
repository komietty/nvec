#ifndef NVEC_FACE_ROSY_FILED_H
#define NVEC_FACE_ROSY_FILED_H
#include "face_tangent_field.h"
#include "linearsolver/linear_solver.h"

namespace pddg {
class FaceRosyField : public FaceTangentField {
public:

    FaceRosyField(const Hmesh& mesh, const int nRosy): FaceTangentField(mesh, nRosy) {}
    FaceRosyField(const Hmesh& mesh, const int nRosy, FieldType type): FaceTangentField(mesh, nRosy) {
        SprsC L = connectionLaplacian();
        SprsC M = galerkinMassMatrix();
        switch (type) {
            case FieldType::Smoothest: {
                compressed = solveSmallestEig(L, M);
                break;
            }
            case FieldType::CurvatureAligned: {
                assert(rosyN == 2 || rosyN == 4);
                constexpr double lambda = 0;
                VecXc D = principalCurvatureDir();
                if (rosyN == 4) D = D.array().square();
                VecXc rhs = M * D / sqrt(abs((D.adjoint() * M * D)[0]));
                SprsC lhs = L - lambda * M;
                compressed = solveSquare(lhs, rhs);
                break;
            }
            default: { break; }
        }
        auto r = std::polar(1., TwoPI / rosyN);
        field.resize(mesh.nF, rosyN);
        for (auto f: mesh.faces) {
        for (int i = 0; i < rosyN; i++) {
            field(f.id, i) = pow(std::polar(1., std::arg(compressed(f.id))), 1. / rosyN) * pow(r, i);
        }}
    }

    void computeSingular(const VecXd& effort) {
        singular = VecXi::Ones(mesh.nV);
        for (auto v: mesh.verts | std::views::filter([](auto _v) { return !_v.isBoundary(); })) {
            double sum = mesh.angleDefect[v.id] * rosyN;
            for (auto h: v.adjHalfs()) sum += (h.isCanonical() ? -1 : 1) * effort[h.edge().id];
            singular[v.id] = static_cast<int>(std::round(sum / TwoPI));
        }
    }

    void computeConnection() {
        connection = VecXd::Zero(mesh.nE);
        for (auto e: mesh.edges | std::views::filter([](auto _e) { return !_e.isBoundary(); })) {
            Vec3d v = e.half().vec().normalized();
            complex ef(v.dot(e.face0().basisX()), v.dot(e.face0().basisY()));
            complex eg(v.dot(e.face1().basisX()), v.dot(e.face1().basisY()));
            connection(e.id) = eg / ef;
        }
    }

    void computeMatching(const MatchingType type) {
        computeConnection();
        matching = VecXi::Constant(mesh.nE, -1);
        VecXd effort = VecXd::Zero(mesh.nE);
        switch (type) {
            case MatchingType::Principal: { pMatching(effort); break; }
            case MatchingType::Curl:      { cMatching(effort); break; }
            default: { break; }
        }
        computeSingular(effort);
    }

private:
    void pMatching(VecXd& effort) {
        for (auto e: mesh.edges | std::views::filter([](auto e_) { return !e_.isBoundary(); })) {
            auto conn = connection[e.id];
            auto minRot = 1000.;
            auto offset = 0;
            complex coef(1, 0);
            complex transport0 = conn * field(e.face0().id, 0);
            for (int i = 0; i < rosyN; i++) {
                complex jf = field(e.face0().id, i);
                complex jg = field(e.face1().id, i);
                coef *= jg / (conn * jf);
                double r = std::arg(jg / transport0);
                if (abs(r) < abs(minRot)) {
                    offset = i;
                    minRot = r;
                }
            }
            effort[e.id] = std::arg(coef);

            double residuals = 0;
            for (int i = 0; i < rosyN; i++) {
                auto v0 = field(e.face0().id, i);
                auto v1 = field(e.face1().id, (i + offset) % rosyN);
                residuals += std::arg(v1 / (conn * v0));
            }
            matching[e.id] = offset - std::round((residuals - effort[e.id]) / TwoPI);
            //std::cout << residuals - effort[e.id] << std::endl;
            //matching[e.id] = offset;
        }
    }

    void cMatching(VecXd& effort) {
        for (auto e: mesh.edges | std::views::filter([](auto e_) { return !e_.isBoundary(); })) {
            int iM = 0;
            double min = 32767000.;
            Vec3d v = e.half().vec().normalized();
            Face f0 = e.face0();
            Face f1 = e.face1();

            for (int i = 0; i < rosyN; i++) {
                double curl = 0;
                for (int j = 0; j < rosyN; j++) {
                    complex c0 = field(f0.id, j);
                    complex c1 = field(f1.id, (i + j) % rosyN);
                    Row3d v0 = c0.real() * f0.basisX() + c0.imag() * f0.basisY();
                    Row3d v1 = c1.real() * f1.basisX() + c1.imag() * f1.basisY();
                    curl += pow(v.dot(v1 - v0), 2.);
                }
                if (curl < min) { iM = i; min = curl; }
            }

            matching[e.id] = iM;

            complex coef(1, 0);
            for (int i = 0; i < rosyN; i++) {
                auto v0 = field(e.face0().id, i);
                auto v1 = field(e.face1().id, (i + iM) % rosyN);
                coef *= v1 / (v0 * connection[e.id]);
            }
            effort[e.id] = arg(coef);
        }
    }
};
}
#endif
