#include "../src/dumbbell_model.h"
#include <iostream>
#include <stdexcept>

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

int main() {
    CONNECTIVITY valid(8, 2);
    require(DUMBBELL::generate_dumbbell_connectivity(valid) == 0, "valid dumbbell rejected");
    for (long i = 0; i < valid.Np; ++i) {
        const long partner = (i + 4) % 8;
        require(valid.HASH[i](0) == i, "self slot changed");
        require(valid.HASH[i](1) == partner, "partner mapping");
        require(valid.HASH[partner](1) == i, "asymmetric pair");
        require(valid.TOKEN[i] == 2, "token must include self and partner");
    }
    CONNECTIVITY too_small(2, 1);
    require(DUMBBELL::generate_dumbbell_connectivity(too_small) == -1, "capacity not rejected");
    require(too_small.TOKEN[0] == 1, "invalid input partially mutated");
    CONNECTIVITY odd(3, 2);
    require(DUMBBELL::generate_dumbbell_connectivity(odd) == -1, "odd Np not rejected");
    std::cout << "PASS: dumbbell pairs and invalid capacity/count guards\n";
}
