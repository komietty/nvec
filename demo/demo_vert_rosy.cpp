#include "nvec/vert/rosy_field.h"
#include <igl/readOBJ.h>
#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
using namespace pddg;

int main(int argc, char *argv[]) {
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    MatXd V;
    MatXi F;
    igl::readOBJ("/Users/saki/dev/models/gargoyle.obj", V, F);

    auto hmsh = std::make_unique<Hmesh>(V, F);

    polyscope::init();
    int rosyN = 4;
    auto surf = polyscope::registerSurfaceMesh("mesh", hmsh->pos, hmsh->idx);
    auto vf = std::make_unique<VertRosyField>(*hmsh, rosyN, FieldType::CurvatureAligned);
    MatXd rosy(hmsh->nV, 2);
    for (Vert v: hmsh->verts) {
        complex c = std::pow(vf->compressed[v.id], 1. / rosyN);
        rosy.row(v.id) = Row2d(c.real(), c.imag()).normalized();
    }
    auto tq = surf->addVertexTangentVectorQuantity("rosyf", rosy, hmsh->vertBasisX, hmsh->vertBasisY, rosyN);
    tq->setEnabled(true);
    tq->setVectorLengthScale(0.004);

    polyscope::show();
}
