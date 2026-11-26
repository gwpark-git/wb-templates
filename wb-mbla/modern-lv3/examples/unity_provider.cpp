#include "case_driver.h"

#include <cmath>
#include <stdexcept>

namespace {

// The constant-law control accepts any finite RK trial, including undershoots.
class UnityTransport final : public mbla::EquilibriumTransport {
public:
    void validate_case(double freeze_fraction, double reference_peclet) const override {
        if (!(freeze_fraction > 0. && freeze_fraction <= 1.) ||
            !std::isfinite(reference_peclet) || reference_peclet < 0.) {
            throw std::invalid_argument("invalid freeze fraction or reference particle Peclet number");
        }
    }

    mbla::TransportValues evaluate(double volume_fraction) override {
        if (!std::isfinite(volume_fraction)) {
            throw std::invalid_argument("unity transport requires a finite trial concentration");
        }
        return {1., 1.};
    }
};

} // namespace

// The application owns the provider for the duration of the existing case driver.
int main(int argc, char* argv[]) {
    UnityTransport transport;
    return mbla::run_cases(argc, argv, &transport);
}
