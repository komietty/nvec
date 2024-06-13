#ifndef NVEC_JSONIO_H
#define NVEC_JSONIO_H

#include <nlohmann/json.hpp>
using json = nlohmann::json;


namespace pddg {
struct person {
    std::string name;
    std::string address;
    int age;
};

void parse() {
    person p = {"Ned Flanders", "744 Evergreen Terrace", 60};
    json j;
    j["name"] = p.name;
    j["address"] = p.address;
    j["age"] = p.age;
}
}


#endif
