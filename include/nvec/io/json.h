#ifndef NVEC_JSONIO_H
#define NVEC_JSONIO_H
#include <nlohmann/json.hpp>
#include <nvec/face_rosy_field.h>

namespace pddg {
void toJson(const std::string& path, const FaceRosyField& f) {
    using json = nlohmann::json;
    json j;
    std::vector<double> vec(f.field.rows() * 3 * f.rosyN);

    MatXd rawExt(f.mesh.nF, 3 * f.rosyN);
    for (Face face: f.mesh.faces) {
        complex rc0 = f.field(face.id, 0);
        complex rc1 = f.field(face.id, 1);
        complex rc2 = f.field(face.id, 2);
        complex rc3 = f.field(face.id, 3);
        Row3d v0 = (rc0.real() * face.basisX() + rc0.imag() * face.basisY()).normalized();
        Row3d v1 = (rc1.real() * face.basisX() + rc1.imag() * face.basisY()).normalized();
        Row3d v2 = (rc2.real() * face.basisX() + rc2.imag() * face.basisY()).normalized();
        Row3d v3 = (rc3.real() * face.basisX() + rc3.imag() * face.basisY()).normalized();
        // rosyN == 4 case
        vec[face.id * 3 * f.rosyN + 0]  = v0.x();
        vec[face.id * 3 * f.rosyN + 1]  = v0.y();
        vec[face.id * 3 * f.rosyN + 2]  = v0.z();
        vec[face.id * 3 * f.rosyN + 3]  = v1.x();
        vec[face.id * 3 * f.rosyN + 4]  = v1.y();
        vec[face.id * 3 * f.rosyN + 5]  = v1.z();
        vec[face.id * 3 * f.rosyN + 6]  = v2.x();
        vec[face.id * 3 * f.rosyN + 7]  = v2.y();
        vec[face.id * 3 * f.rosyN + 8]  = v2.z();
        vec[face.id * 3 * f.rosyN + 9]  = v3.x();
        vec[face.id * 3 * f.rosyN + 10] = v3.y();
        vec[face.id * 3 * f.rosyN + 11] = v3.z();
    }
    j["BaseType"] = "face";
    j["PolyType"] = "rosy";
    j["NumPoly"] = f.rosyN;
    j["NumFace"] = f.mesh.nF;
    j["Elements"] = vec;

    std::ofstream o(path.c_str());
    o << j.dump() << std::endl;
    o.close();
}

MatXd fromJson(const std::string& path) {
    using json = nlohmann::json;
    using MatXd_RM = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    if (std::ifstream ifs(path); ifs.good()) {
        json j;
        ifs >> j;
        const int nP = j["NumPoly"];
        const int nF = j["NumFace"];
        std::vector<double> elm = j["Elements"];
        assert(elm.size() == nF * nP * 3);
        return Eigen::Map<MatXd_RM>(elm.data(), nF, nP * 3);
    }
}
}
#endif
