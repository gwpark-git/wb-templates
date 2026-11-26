#pragma once

namespace mbla {

struct TransportValues {
    double collective_diffusivity_over_D0;
    double viscosity_over_solvent;
};

// Runtime law in the local-equilibrium limit. The caller owns its lifetime.
class EquilibriumTransport {
public:
    virtual ~EquilibriumTransport() = default;
    virtual void validate_case(double freeze_fraction, double reference_peclet) const = 0;
    virtual TransportValues evaluate(double volume_fraction) = 0;
};

// Standard CLI dispatch. Only the paired LIVE selections use an injected law.
int run_cases(int argc, char* argv[], EquilibriumTransport* transport = nullptr);

} // namespace mbla
