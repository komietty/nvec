#ifndef NVEC_FACE_COMBING_H
#define NVEC_FACE_COMBING_H
#include <ranges>
#include "face_rosy_field.h"

namespace pddg {
inline std::unique_ptr<FaceRosyField> computeComb(
        const FaceRosyField& field,
        const std::vector<bool>& seam
) {

    auto N = field.rosyN;
    auto F = std::make_unique<FaceRosyField>(field.mesh, field.rosyN);
    F->field.resize(field.mesh.nF, N);
    F->connection = field.connection;
    F->compressed = field.compressed;
    F->singular = field.singular;
    VecXi turns = VecXi::Zero(F->mesh.nF);
    VecXi visit = VecXi::Zero(F->mesh.nF);
    std::queue<std::pair<int, int>> matchingQ;
    matchingQ.emplace(0, 0);

    do {
        std::pair<int, int> pop = matchingQ.front();
        matchingQ.pop();
        auto [i, m] = pop;
        int r = N - m;
        if (visit(i) == 1) continue;
        visit(i) = 1;
        turns(i) = m;
        F->field.block(i, 0, 1, r) = field.field.block(i, m, 1, r);
        F->field.block(i, r, 1, m) = field.field.block(i, 0, 1, m);
        Face face = F->mesh.faces[i];

        for (Half h: face.adjHalfs()) {
            Edge e = h.edge();
            bool b = e.face0().id == i;
            Face fn = b ? e.face1() : e.face0();
            int n = (b ? 1 : -1) * field.matching[e.id];
            n = (n + m + 1000 * N) % N;
            if (!e.isBoundary() && !visit(fn.id) && !seam[e.id]) matchingQ.emplace(fn.id, n);
        }
    } while (!matchingQ.empty());

    F->matching = VecXi::Constant(F->mesh.nE, -1);
    for (Edge e: F->mesh.edges | std::views::filter([](auto e) { return !e.isBoundary();}))
        F->matching[e.id] = (turns[e.face0().id] - turns[e.face1().id] + field.matching[e.id] + 1000 * N) % N;

    return F;
}
}

#endif
