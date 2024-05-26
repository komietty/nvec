#ifndef NVEC_VERT_ROSY_FIELD_H
#define NVEC_VERT_ROSY_FIELD_H
#include "base_field.h"
#include "linearsolver/linear_solver.h"

namespace pddg {
class VertRosyField : public BaseVectorField {
public:
    VertRosyField(const Hmsh& m, const int n, const FieldType type): BaseVectorField(m, n) {
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
                VecXc D = principalCurvatureDirection();
                if (rosyN == 4) D = D.array().square();
                VecXc rhs = M * D / sqrt(abs((D.adjoint() * M * D)[0]));
                SprsC lhs = L - lambda * M;
                compressed = solveSquare(lhs, rhs);
                break;
            }
            default: { break; }
        }
    }

    SprsC connectionLaplacian() const {
        SprsC S(mesh.nV, mesh.nV);
        SprsC I(mesh.nV, mesh.nV);
        VecXc transport(mesh.nH);
        std::vector<TripC> T;

        for (Edge e: mesh.edges) {
            Half h1 = e.half();
            Half h2 = e.half().twin();
            auto v1 = std::polar(1., h1.varg());
            auto v2 = std::polar(1., h2.varg());
            auto r = -v2 / v1;
            transport[h1.id] = r;
            transport[h2.id] = complex(1, 0) / r;
        }

        for (Half h: mesh.halfs) {
            auto w = h.edge().cot();
            auto r = transport[h.twin().id];
            T.emplace_back(h.tail().id, h.tail().id, w);
            T.emplace_back(h.tail().id, h.head().id, -w * pow(r, rosyN));
        }

        S.setFromTriplets(T.begin(), T.end());
        I.setIdentity();
        return S + 1e-9 * I;
    }

    SprsC galerkinMassMatrix() const {
        SprsC S(mesh.nV, mesh.nV);
        std::vector<TripC> T;
        for (Face f: mesh.faces) {
            auto a = f.area();
            auto h = f.half();
            std::array vs{h.tail(), h.next().tail(), h.prev().tail()};
            for (int i = 0; i < 3; i++) {
                int iv = vs[i].id;
                int jv = vs[(i + 1) % 3].id;
                int kv = vs[(i + 2) % 3].id;
                T.emplace_back(iv, iv, complex{a / 6., 0});
                T.emplace_back(iv, jv, complex{a / 12., 0});
                T.emplace_back(iv, kv, complex{a / 12., 0});
            }
        }
        S.setFromTriplets(T.begin(), T.end());
        return S;
    }

    VecXc principalCurvatureDirection() const {
        VecXc D(mesh.nV);
        for (Vert v: mesh.verts) {
            complex dir{0, 0};
            for (Half h: v.adjHalfs()) {
                auto l = h.len();
                auto c = std::polar(1., h.varg()) * l;
                dir += -c * c / l * h.darg();
            }
            D[v.id] = dir * 0.25;
        }
        return D;
    }

    VecXi computeSingularNum();

    static MatXd convert_to_extrinsic_field(const Mesh& m, const VecXc& f, const int nSym) {
        MatXd ext(m.nV, 3 * nSym);
        for (Vert v: m.verts) {
            complex c = std::pow(f[v.id], 0.25);
            complex r = std::exp(complex(0, 1) * TwoPI / (double) nSym);
            for (int j = 0; j < nSym; j++) {
                ext.block(v.id, j * 3, 1, 3) = (Row3d) (c.real() * v.basisX() + c.imag() * v.basisY());
                c *= r;
            }
        }
        return ext;
    }
};

inline VecXi VertRosyField::computeSingularNum() {
    auto fmodPI = [](double t) { return t - TwoPI * floor((t + PI) / TwoPI); };
    VecXi N(mesh.nF);
    for (Face f: mesh.faces) {
        double omega = 0;
        double idx = 0;
        for (Half h: f.adjHalfs()) {
            double phiI = std::arg(compressed[h.tail().id]);
            double phiJ = std::arg(compressed[h.head().id]);
            double thetaI = h.varg();
            double thetaJ = h.twin().varg() + PI;
            double dTheta = thetaI - thetaJ;
            omega += dTheta;
            idx += fmodPI(phiJ - phiI + rosyN * dTheta);
        }
        idx -= rosyN * fmodPI(omega);
        N[f.id] = std::lround(idx / TwoPI);
    }
    return N;
}
}

#endif
