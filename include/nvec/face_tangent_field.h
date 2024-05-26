#ifndef NVEC_FACE_TANGENT_FIELD_H
#define NVEC_FACE_TANGENT_FIELD_H
#include "base_field.h"

namespace pddg {
class FaceTangentField: public BaseVectorField {
public:
    FaceTangentField(const Hmsh& m, const int n) : BaseVectorField(m, n) { }

    [[nodiscard]] SprsC connectionLaplacian() const {
        SprsC S(mesh.nF, mesh.nF);
        SprsC I(mesh.nF, mesh.nF);
        VecXc transport(mesh.nH);
        std::vector<TripC> T;

        for (Edge e: mesh.edges) {
            Half h1 = e.half();
            Half h2 = e.half().twin();
            auto v1 = std::polar(1., h1.farg());
            auto v2 = std::polar(1., h2.farg());
            auto r = -v2 / v1;
            transport[h1.id] = r;
            transport[h2.id] = complex(1, 0) / r;
        }

        for (Face f: mesh.faces) {
            for (Half h: f.adjHalfs()) {
                if (h.twin().isBoundary()) continue;
                auto w = 1.;
                auto r = transport[h.twin().id];
                Face twinF = h.twin().face();
                T.emplace_back(f.id, f.id, w);
                T.emplace_back(f.id, twinF.id, -w * pow(r, rosyN));
            }}

        S.setFromTriplets(T.begin(), T.end());
        I.setIdentity();
        return S + 1e-9 * I;
    }

    [[nodiscard]] SprsC galerkinMassMatrix() const {
        SprsC S(mesh.nF, mesh.nF);
        std::vector<TripC> T;
        for (Face f: mesh.faces) { T.emplace_back(f.id, f.id, f.area()); }
        S.setFromTriplets(T.begin(), T.end());
        return S;
    }

    [[nodiscard]] VecXc principalCurvatureDir() const {
        VecXc D(mesh.nF);
        for (Face f: mesh.faces) {
            complex dir{0, 0};
            for (Half h: f.adjHalfs()) {
                auto l = h.len();
                auto c = std::polar(1., h.farg()) * l;
                dir += -c * c / l * h.darg();
            }
            D[f.id] = dir * 0.25;
        }
        return D;
    }

};
}

#endif
