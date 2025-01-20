#include <memory>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/frame_field_deformer.h>
#include "hmesh/hmesh.h"
#include "nvec/face/frame_field.h"

/// --------------------------------------------------- ///
/// This demo is not working yet.
/// --------------------------------------------------- ///

using namespace pddg;
MatXd V;
MatXi F;

int main(int argc, char *argv[]) {
    polyscope::options::autocenterStructures = true;

    //igl::readOBJ(argv[1], V, F);
    igl::readOFF(argv[1], V, F);
    auto hmsh = std::make_unique<Hmesh>(V, F);

    MatXd HN;
    SprsD L, M, Minv;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    // Laplace-Beltrami of position
    HN = -Minv * (L * V);
    // Extract magnitude as mean curvature
    VecXd H = HN.rowwise().norm();

    // Compute curvature directions via quadric fitting
    MatXd PD1, PD2;
    VecXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    // mean curvature
    H = 0.5 * (PV1 + PV2);
    VecXd K = PV1.array() * PV2.array();
    MatXd anisoPD1 = PV1.cwiseSqrt().cwiseInverse().asDiagonal() * PD1;
    MatXd anisoPD2 = PV2.cwiseSqrt().cwiseInverse().asDiagonal() * PD2;

    //--- frame field bgn ---//
    VecXi b(hmsh->nF);
    MatXd bc1(hmsh->nF, 3);
    MatXd bc2(hmsh->nF, 3);
    MatXd ext(hmsh->nF, 3);
    MatXd FF1, FF2;
    MatXd V_deformed;
    MatXd FF1_deformed;
    MatXd FF2_deformed;

    for (Face f: hmsh->faces) {
        b[f.id] = f.id;
        int id0 = f.half().tail().id;
        int id1 = f.half().next().tail().id;
        int id2 = f.half().prev().tail().id;

        Row3d v1 = anisoPD1.row(id0);
        Row3d v2 = anisoPD2.row(id0);
        Row3d v1n = anisoPD1.row(id1);
        Row3d v2n = anisoPD2.row(id1);
        Row3d v1p = anisoPD1.row(id2);
        Row3d v2p = anisoPD2.row(id2);

        if (v1n.dot(v1) < 0) v1n *= -1;
        if (v2n.dot(v2) < 0) v2n *= -1;
        if (v1p.dot(v1) < 0) v1p *= -1;
        if (v2p.dot(v2) < 0) v2p *= -1;

        bc1.row(f.id) = (v1 + v1n + v1p) / 3.;
        bc2.row(f.id) = (v2 + v2n + v2p) / 3.;
        ext.row(f.id) = (v1 + v1n + v1p).normalized();
    }

    frame_field(V, F, b, bc1, bc2, ext, FF1, FF2);
    igl::frame_field_deformer(V, F, FF1, FF2, V_deformed, FF1_deformed, FF2_deformed);

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    auto surf = polyscope::registerSurfaceMesh("mesh", hmsh->pos, hmsh->idx);

    /*--- visuailize curvature ---*/
    surf->setSurfaceColor({0, 10./ 255., 27./ 255.});
    surf->addVertexScalarQuantity("Mean Curvature libigl", H);
    surf->addVertexScalarQuantity("Gaussian Curvature libigl", K);
    surf->addVertexScalarQuantity("Principal Curvature libigl Max", PV1);
    surf->addVertexScalarQuantity("Principal Curvature libigl Min", PV2);
    surf->addVertexVectorQuantity("Principal Curvature Max", anisoPD1);
    surf->addVertexVectorQuantity("Principal Curvature Min", anisoPD2);
    surf->addFaceVectorQuantity("bc1", bc1);
    surf->addFaceVectorQuantity("bc2", bc2);
    surf->resetTransform();
    surf->setSmoothShade(true);

    auto defo = polyscope::registerSurfaceMesh("defo", V_deformed, hmsh->idx);
    auto hmsh2 = std::make_unique<Hmesh>(V_deformed, F);
    defo->addFaceVectorQuantity("defo_FF1", FF1_deformed);
    defo->addFaceVectorQuantity("defo_FF2", FF2_deformed);
    defo->addVertexScalarQuantity("Gaussian Curvature", hmsh2->angleDefect);

    polyscope::show();
}
