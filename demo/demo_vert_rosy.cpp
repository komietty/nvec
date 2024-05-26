
int main(int argc, char *argv[]) {
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    MatXd V;
    MatXi F;
    igl::readOBJ("/Users/komietty/dev/models/fandisk.obj", V, F);
    auto mesh = std::make_unique<Mesh>(V, F);
    auto tree = std::make_unique<Tree>(*mesh);
    tree->calcHomologyGens(false);

    polyscope::init();
    int rosyN = 4;
    auto surf = polyscope::registerSurfaceMesh("mesh", mesh->pos, mesh->idx);
    auto vf = std::make_unique<VertRosyField>(*mesh, rosyN, FieldType::CurvatureAligned);
    //ff->computeMatching(MatchingType::Principal);
    MatXd rosy(mesh->nV, 2);
    for (Vert v: mesh->verts) {
        complex c = std::pow(vf->compressed[v.id], 1. / rosyN);
        rosy.row(v.id) = Row2d(c.real(), c.imag()).normalized();
    }
    auto tq = surf->addVertexTangentVectorQuantity("rosyf", rosy, mesh->vertBasisX, mesh->vertBasisY, rosyN);
    tq->setEnabled(true);
    tq->setVectorLengthScale(0.004);

    //std::vector<glm::vec3> singPos;
    //std::vector<double> singVal;
    //for (Face f: mesh->faces) {
    //    if (int s = vf->singular[f.id]; s != 0) {
    //        Row3d p = f.center();
    //        singPos.emplace_back(p.x(), p.y(), p.z());
    //        singVal.emplace_back(s);
    //    }
    //}
    //auto pc = polyscope::registerPointCloud("face field singulars", singPos);
    //pc->addScalarQuantity("face field singular nums", singVal);
    //pc->setEnabled(true);
    //pc->setPointRadius(0.005);
    //pc->resetTransform();

    polyscope::show();
}
