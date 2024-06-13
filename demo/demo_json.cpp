#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;  // 推奨されているエイリアス

namespace ns {
struct person {
    std::string name;
    std::string address;
    int age;
};
}

int main() {
    ns::person p = {"Ned Flanders", "744 Evergreen Terrace", 60};

    std::vector<int> c_vector;
    c_vector.resize(1000000);
    for(int i = 0; i < c_vector.size(); i++) {
        c_vector[i] = i;
    }
    json j_vec(c_vector);

    json j;
    j["name"] = p.name;
    j["address"] = p.address;
    j["age"] = p.age;
    j["list"] = j_vec;

    std::ofstream o("pretty.json");
    o << std::setw(4) << j << std::endl;
}
