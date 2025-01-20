#ifndef NVEC_BASE_FIELD_H
#define NVEC_BASE_FIELD_H
#include "hmesh/hmesh.h"
#include <ranges>

namespace pddg {
class BaseVectorField {
public:
    const Hmesh& mesh;
    const int rosyN;
    MatXc field;
    VecXc connection;
    VecXc compressed;
    VecXi matching;
    VecXi singular;
    BaseVectorField(const Hmesh& m, const int n) : mesh(m), rosyN(n) { }
};

enum class FieldType { UnSpecified, Smoothest, CurvatureAligned };
enum class MatchingType { Principal, Curl };
}

#endif
