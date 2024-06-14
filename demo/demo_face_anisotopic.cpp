#include "nvec/Face/FrameField.h"
#include "nvec/face_rosy_field.h"
#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/frame_field_deformer.h>
#include <igl/frame_to_cross_field.h>
#include <igl/readDMAT.h>
#include <igl/rotate_vectors.h>

using namespace pddg;

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;

// Input frame field constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc1;
Eigen::MatrixXd bc2;

// Interpolated frame field
Eigen::MatrixXd FF1, FF2;

// Deformed mesh
Eigen::MatrixXd V_deformed;
Eigen::MatrixXd B_deformed;

// Frame field on deformed
Eigen::MatrixXd FF1_deformed;
Eigen::MatrixXd FF2_deformed;

// Cross field on deformed
Eigen::MatrixXd X1_deformed;
Eigen::MatrixXd X2_deformed;

/*
bool key_down(
        igl::opengl::glfw::Viewer& viewer,
        unsigned char key,
        int modifier
        )
{
    using namespace std;
    using namespace Eigen;

    if (key <'1' || key >'6')
        return false;

    viewer.data().clear();
    viewer.data().show_lines = false;
    viewer.data().show_texture = false;

    if (key == '1')
    {
        // Frame field constraints
        viewer.data().set_mesh(V, F);

        MatrixXd F1_t = MatrixXd::Zero(FF1.rows(),FF1.cols());
        MatrixXd F2_t = MatrixXd::Zero(FF2.rows(),FF2.cols());
        // Highlight in red the constrained faces
        MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
        for (unsigned i=0; i<b.size();++i)
        {
            C.row(b(i)) << 1, 0, 0;
            F1_t.row(b(i)) = bc1.row(i);
            F2_t.row(b(i)) = bc2.row(i);
        }

        viewer.data().set_colors(C);

        MatrixXd C1,C2;
        VectorXd K1 = F1_t.rowwise().norm();
        VectorXd K2 = F2_t.rowwise().norm();
        igl::jet(K1,true,C1);
        igl::jet(K2,true,C2);

        viewer.data().add_edges(B - global_scale*F1_t, B + global_scale*F1_t ,C1);
        viewer.data().add_edges(B - global_scale*F2_t, B + global_scale*F2_t ,C2);
    }

    if (key == '2')
    {
        // Frame field
        viewer.data().set_mesh(V, F);
        MatrixXd C1,C2;
        VectorXd K1 = FF1.rowwise().norm();
        VectorXd K2 = FF2.rowwise().norm();
        igl::jet(K1,true,C1);
        igl::jet(K2,true,C2);

        viewer.data().add_edges(B - global_scale*FF1, B + global_scale*FF1 ,C1);
        viewer.data().add_edges(B - global_scale*FF2, B + global_scale*FF2 ,C2);

        // Highlight in red the constrained faces
        MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
        for (unsigned i=0; i<b.size();++i)
            C.row(b(i)) << 1, 0, 0;
        viewer.data().set_colors(C);

    }

    if (key == '3')
    {
        // Deformed with frame field
        viewer.data().set_mesh(V_deformed, F);
        viewer.data().add_edges(B_deformed - global_scale*FF1_deformed, B_deformed + global_scale*FF1_deformed ,Eigen::RowVector3d(1,0,0));
        viewer.data().add_edges(B_deformed - global_scale*FF2_deformed, B_deformed + global_scale*FF2_deformed ,Eigen::RowVector3d(0,0,1));
        viewer.data().set_colors(RowVector3d(1,1,1));
    }

    if (key == '4')
    {
        // Deformed with cross field
        viewer.data().set_mesh(V_deformed, F);
        viewer.data().add_edges(B_deformed - global_scale*X1_deformed, B_deformed + global_scale*X1_deformed ,Eigen::RowVector3d(0,0,1));
        viewer.data().add_edges(B_deformed - global_scale*X2_deformed, B_deformed + global_scale*X2_deformed ,Eigen::RowVector3d(0,0,1));
        viewer.data().set_colors(RowVector3d(1,1,1));
    }

    // Replace the standard texture with an integer shift invariant texture
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;
    line_texture(texture_R, texture_G, texture_B);
    viewer.data().set_texture(texture_R, texture_B, texture_G);
    viewer.core().align_camera_center(viewer.data().V,viewer.data().F);

    return false;
}
*/

int main(int argc, char *argv[]) {
    // Load a mesh in OBJ format
    igl::readOBJ("/Users/komietty/dev/models/bumpy-cube.obj", V, F);
    igl::barycenter(V, F, B);

    // Compute scale for visualizing fields
    global_scale = .2 * igl::avg_edge_length(V, F);

    // Load constraints
    MatXd temp;
    igl::readDMAT("/Users/komietty/dev/models/bumpy-cube.dmat",temp);

    b   = temp.block(0,0,temp.rows(),1).cast<int>();
    bc1 = temp.block(0,1,temp.rows(),3);
    bc2 = temp.block(0,4,temp.rows(),3);

    /// ------
    auto mesh = std::make_unique<Hmsh>(V, F);

    int rosyN = 4;
    auto rawf = std::make_unique<FaceRosyField>(*mesh, rosyN, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    MatXd rawExt(mesh->nF, 3);
    for (Face f: mesh->faces) {
        complex c1 = rawf->field(f.id, 0);
        //complex c1 = rawf->compressed(f.id);
        rawExt.row(f.id) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
    }
    /// -----

    // Interpolate the frame field
    pddg::frame_field(V, F, b, bc1, bc2, rawExt, FF1, FF2);
    std::cout << "1" << std::endl;


    // Deform the mesh to transform the frame field in a cross field
    igl::frame_field_deformer(V,F,FF1,FF2,V_deformed,FF1_deformed,FF2_deformed);
    std::cout << "2" << std::endl;

    polyscope::init();
    auto msh2 = std::make_unique<Hmsh>(V_deformed, F);
    auto surf = polyscope::registerSurfaceMesh("mesh", msh2->pos, msh2->idx);
    surf->setSurfaceColor({0, 10./ 255., 27./ 255.});
    surf->setSmoothShade(true);
    auto f1 = surf->addFaceVectorQuantity("raw ext 1", FF1_deformed);
    auto f2 = surf->addFaceVectorQuantity("raw ext 2", FF2_deformed);
    auto f3 = surf->addFaceVectorQuantity("raw ext 3", -FF1_deformed);
    auto f4 = surf->addFaceVectorQuantity("raw ext 4", -FF2_deformed);
    f1->setVectorLengthScale(0.004)->setEnabled(true);
    f2->setVectorLengthScale(0.004)->setEnabled(true);
    f3->setVectorLengthScale(0.004)->setEnabled(true);
    f4->setVectorLengthScale(0.004)->setEnabled(true);
    polyscope::show();

    // Compute face barycenters deformed mesh
    igl::barycenter(V_deformed, F, B_deformed);
    std::cout << "3" << std::endl;

    // Find the closest crossfield to the deformed frame field
    igl::frame_to_cross_field(V_deformed,F,FF1_deformed,FF2_deformed, X1_deformed);
    std::cout << "4" << std::endl;

    // Find a smooth crossfield that interpolates the deformed constraints
    MatXd bc_x(b.size(),3);
    for (unsigned i=0; i<b.size();++i)
        bc_x.row(i) = X1_deformed.row(b(i));
}