#include "hs_mc.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
using State = std::unique_ptr<HS, decltype(&hs_destroy)>;

void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
}

int main() {
    try {
        const double pair[] = {0., 0., 0., 1., 0., 0.};
        const State touching(hs_create(2, 5, 42, 54, 0.15, 0, pair), hs_destroy);
        check(bool(touching), "exact contact should be allowed");
        const double across_boundary[] = {4.5, 0., 0.};
        check(!hs_allowed(touching.get(), 1, across_boundary), "PBC overlap not rejected");
        const double overlapping[] = {0., 0., 0., 4.5, 0., 0.};
        const State invalid(hs_create(2, 5, 42, 54, 0.15, 0, overlapping), hs_destroy);
        check(!invalid, "initial PBC overlap not rejected");

        const double outside_contact[] = {0., 0., 0., 1.015, 0., 0.};
        const State shell(hs_create(2, 5, 42, 54, 0.15, 0, outside_contact), hs_destroy);
        check(bool(shell), "pair fixture initialization");
        double measured[7];
        hs_measure(shell.get(), measured);
        check(std::abs(measured[0] - 1.015 * 1.015) < 1e-14, "wrong squared separation");
        check(measured[1] == 0 && measured[2] == 0 && measured[3] == 0 &&
              measured[4] == 1 && measured[5] == 1, "cumulative unordered shell counts");
        const State ideal(hs_create(1, 5, 42, 54, 0.15, 0, pair), hs_destroy);
        check(bool(ideal) && hs_sweeps(ideal.get(), 30) == 30, "single-particle acceptance");
        hs_measure(ideal.get(), measured);
        check(measured[0] == 1.05 * 1.05 && measured[5] == 0, "ideal capped separation");

        std::vector<double> initial;
        for (double x : {0.5, 1.9, 3.3})
            for (double y : {0.5, 1.9, 3.3})
                for (double z : {0.5, 1.9, 3.3})
                    initial.insert(initial.end(), {x, y, z});
        const State grid(hs_create(27, 4.2, 91371, 11, 0.3, 0, initial.data()), hs_destroy);
        const State brute(hs_create(27, 4.2, 91371, 11, 0.3, 1, initial.data()), hs_destroy);
        check(bool(grid) && bool(brute), "trajectory fixture initialization");
        std::array<double, 81> a, b, da, db;
        unsigned accepted_total = 0;
        for (unsigned iteration = 0; iteration < 100; ++iteration) {
            const auto accepted = hs_sweeps(grid.get(), 5);
            check(accepted == hs_sweeps(brute.get(), 5), "cell/brute acceptance correlation");
            accepted_total += accepted;
            hs_copy(grid.get(), a.data(), da.data());
            hs_copy(brute.get(), b.data(), db.data());
            check(a == b && da == db && hs_check_lists(grid.get()), "cell/brute trajectory or grid drift");
        }
        check(accepted_total > 0 && accepted_total < 100 * 5 * 27, "need accepted and rejected moves");
        std::cout << "PASS contact, periodic overlap, shell counts, ideal motion and paired linked-cell trajectory\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
