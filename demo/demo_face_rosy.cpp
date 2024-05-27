#include "nvec/face_combing.h"
#include "nvec/face_seaming.h"
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
    igl::readOBJ("/Users/komietty/dev/models/gargoyle.obj", V, F);
    auto mesh = std::make_unique<Hmsh>(V, F);

    polyscope::init();
    int rosyN = 4;
    auto surf = polyscope::registerSurfaceMesh("mesh", mesh->pos, mesh->idx);
    surf->setSurfaceColor({0, 10./ 255., 27./ 255.});
    surf->setSmoothShade(true);
    auto rawf = std::make_unique<FaceRosyField>(*mesh, rosyN, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = computeSeam(*rawf);
    auto cmbf = computeComb(*rawf, seam);
    MatXd rawInt(mesh->nF, 2);
    MatXd cmbInt(mesh->nF, 2);
    MatXd rawExt(mesh->nF, 3);
    MatXd cmbExt(mesh->nF, 3);
    for (Face f: mesh->faces) {
        complex c1 = rawf->field(f.id, 0);
        complex c2 = cmbf->field(f.id, 0);
        rawInt.row(f.id) = Row2d(c1.real(), c1.imag()).normalized();
        cmbInt.row(f.id) = Row2d(c2.real(), c2.imag()).normalized();
        rawExt.row(f.id) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
        cmbExt.row(f.id) = (c2.real() * f.basisX() + c2.imag() * f.basisY()).normalized();
    }
    auto rawFQ = surf->addFaceVectorQuantity("raw ext", rawExt);
    auto cmbFQ = surf->addFaceVectorQuantity("cmb ext", cmbExt);
    rawFQ->setEnabled(false);
    cmbFQ->setEnabled(true);
    rawFQ->setVectorLengthScale(0.004);
    cmbFQ->setVectorLengthScale(0.004);

    std::vector<glm::vec3> singPos;
    std::vector<double> singVal;
    for (Vert v: mesh->verts) {
        if (int s = rawf->singular[v.id]; s != 0) {
            Row3d p = v.pos();
            singPos.emplace_back(p.x(), p.y(), p.z());
            singVal.emplace_back(s);
        }
    }
    auto pc = polyscope::registerPointCloud("face field singulars", singPos);
    pc->addScalarQuantity("face field singular nums", singVal);
    pc->setEnabled(true);
    pc->setPointRadius(0.005);
    pc->resetTransform();

    /*--- visuailize seam ---*/
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    size_t counter = 0;
    for (auto e: mesh->edges) {
        if(seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            nodes.emplace_back(p1.x(), p1.y(), p1.z());
            nodes.emplace_back(p2.x(), p2.y(), p2.z());
            edges.emplace_back(std::array{ counter, counter + 1 });
            counter += 2;
        }
    }
    auto pn = polyscope::registerCurveNetwork("seam", nodes, edges);
    pn->resetTransform();
    pn->setRadius(0.001);

    polyscope::show();
}
