#ifndef NVEC_BASE_FIELD_H
#define NVEC_BASE_FIELD_H
#include "hmsh.h"
#include <ranges>

namespace pddg {
class BaseVectorField {
public:
    const Hmsh& mesh;
    const int rosyN;
    MatXc field;
    VecXc connection;
    VecXc compressed;
    VecXi matching;
    VecXi singular;
    BaseVectorField(const Hmsh& m, const int n) : mesh(m), rosyN(n) { }
};

enum class FieldType { UnSpecified, Smoothest, CurvatureAligned };
enum class MatchingType { Principal, Curl };
}

#endif
