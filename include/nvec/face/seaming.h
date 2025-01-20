#ifndef NVEC_FACE_SEAM_H
#define NVEC_FACE_SEAM_H
#include <igl/dijkstra.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/cut_mesh_from_singularities.h>
#include <set>
#include "hmesh/typedefs.h"
#include "rosy_field.h"

namespace pddg {
inline void cut_mesh_with_singularities(
    const MatXd& V,
    const MatXi& F,
    const std::vector<std::vector<int>>& VF,
    const std::vector<std::vector<int>>& VV,
    const MatXi& TT,
    const MatXi& TTi,
    const VecXi& singularities,
    Eigen::MatrixXi& cuts
) {
    //first, get a spanning tree for the mesh (no missmatch needed)
    igl::cut_mesh_from_singularities(V, F, MatXd::Zero(F.rows(), 3).eval(), cuts);

    std::set<int> vertices_in_cut;
    for (int i = 0; i < cuts.rows(); ++i)
    for (int j = 0; j < cuts.cols(); ++j)
        if (cuts(i, j)) vertices_in_cut.insert(F(i, j));

    //then, add all singularities one by one by using Dijkstra's algorithm
    for (int i = 0; i < singularities.rows(); ++i) {
        std::vector<int> path;
        VecXd min_distance;
        VecXi previous;
        int vertex_found = igl::dijkstra(singularities[i], vertices_in_cut, VV, min_distance, previous);
        if (vertex_found == -1)
            // this means that there are no cuts
            path.push_back(singularities[i]);
        else
            igl::dijkstra(vertex_found, previous, path);

        vertices_in_cut.insert(path.begin(), path.end());

        //insert to cut
        for (int ii = 0; ii < path.size() - 1; ++ii) {
            const int& v0 = path[ii];
            const int& v1 = path[ii + 1];

            std::vector<int> vf0 = VF[v0]; std::sort(vf0.begin(), vf0.end());
            std::vector<int> vf1 = VF[v1]; std::sort(vf1.begin(), vf1.end());

            std::vector<int> common_face_v(std::max(vf0.size(), vf1.size()));
            std::vector<int>::iterator it;
            it = std::set_intersection(vf0.begin(), vf0.end(), vf1.begin(), vf1.end(), common_face_v.begin());
            common_face_v.resize(it - common_face_v.begin());
            assert(common_face_v.size() == 2);

            const int& fi = common_face_v[0];
            int j = -1;
            for (unsigned z = 0; z < 3; ++z)
                if ((F(fi, z) == v0 && F(fi, (z + 1) % 3) == v1) ||
                    (F(fi, z) == v1 && F(fi, (z + 1) % 3) == v0)) { j = z; }
            assert(j != -1);
            cuts(fi, j) = 1;
            cuts(TT(fi, j), TTi(fi, j)) = 1;
        }
    }
}

//Wrapper of the above with only vertices and faces as mesh input
inline void cut_mesh_with_singularities(
    const MatXd& V,
    const MatXi& F,
    const VecXi& singularities,
    MatXi& cuts
) {
    std::vector<std::vector<int>> VF, VFi;
    std::vector<std::vector<int>> VV;
    igl::vertex_triangle_adjacency(V, F, VF, VFi);
    igl::adjacency_list(F, VV);
    MatXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities, cuts);
}

inline std::vector<bool> computeSeam(const FaceRosyField& f) {
    MatXi face2cut;
    std::vector seam(f.mesh.nE, false);
    std::vector<int> svids;
    for (Vert v: f.mesh.verts)
        if (f.singular[v.id] != 0)
            svids.emplace_back(v.id);
    VecXi s = Eigen::Map<VecXi>(svids.data(), svids.size());
    cut_mesh_with_singularities(f.mesh.pos, f.mesh.idx, s, face2cut);
    for (int iF = 0; iF < f.mesh.nF; iF++) {
    for (int j = 0; j < 3; j++) {
        if (face2cut(iF, j)) seam[f.mesh.face2edge(iF, j)] = true;
    }}
    return seam;
}
}
#endif
