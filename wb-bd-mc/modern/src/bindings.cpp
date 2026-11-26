#include "bdmc/core.hpp"
#include "bdmc/spatial.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <string>

namespace py = pybind11;

#ifdef BDMC_HAS_EXTRA_BINDINGS
void bind_extra(py::module_& module);
#endif

namespace {

py::ssize_t array_size(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<py::ssize_t>::max())) {
        throw std::overflow_error("array dimension exceeds Python ssize_t");
    }
    return static_cast<py::ssize_t>(value);
}

template <typename T>
py::array_t<T> copied_array(const std::vector<T>& values,
                           const std::vector<py::ssize_t>& shape) {
    // Owned copies make Python analysis independent of the next native step and
    // of the simulation object's lifetime. A writable NumPy result must never
    // become an undocumented route to mutate the native particle state.
    py::array_t<T> array(shape);
    if (!values.empty()) std::copy(values.begin(), values.end(), array.mutable_data());
    return array;
}

py::array_t<std::uint64_t> bond_array(const std::vector<bdmc::Bond>& bonds) {
    py::array_t<std::uint64_t> result({array_size(bonds.size()), py::ssize_t{2}});
    auto* output = result.mutable_data();
    for (const auto& bond : bonds) {
        *output++ = bond.first;
        *output++ = bond.second;
    }
    return result;
}

py::dict spec_dict(const bdmc::SimulationSpec& spec) {
    py::dict result;
    result["model"] = bdmc::model_name(spec.model);
    result["n_particles"] = spec.n_particles;
    result["dimension"] = spec.dimension;
    result["box_lengths"] = spec.box_lengths;
    result["dt"] = spec.dt;
    result["mobility"] = spec.mobility;
    result["diffusivity"] = spec.diffusivity;
    result["spring_alpha"] = spec.spring_alpha;
    result["seed"] = spec.seed;
    result["noise_policy"] = bdmc::noise_name(spec.noise_policy);
    return result;
}

py::dict metadata(const bdmc::SimulationSpec& spec, bool explicit_initial_positions) {
    auto result = spec_dict(spec);
    result["schema"] = "bdmc.observation.v1";
    result["backend"] = "modern_cpu_reference";
    result["execution"] = "serial CPU; no GPU/MPI backend";
    result["boundary"] = "orthorhombic periodic; minimum image in [-L/2,L/2)";
    result["time_unit"] = "tau_B (dimensionless Brownian time)";
    result["length_unit"] = "dimensionless; caller must preserve the selected legacy normalization";
    result["force_convention"] = "dx_drift=mobility*F*dt; F=-gradient(U)";
    result["noise_increment_variance"] = "2*diffusivity*dt per coordinate";
    result["noise_finite_step_law"] = spec.noise_policy == bdmc::NoisePolicy::LegacyUniform
        ? "uniform[-sqrt(6*diffusivity*dt),+sqrt(6*diffusivity*dt))"
        : "Gaussian(mean=0,variance=2*diffusivity*dt)";
    result["rng"] = "std::mt19937_64; explicit top-53-bit uniform; Box-Muller Gaussian";
    result["rng_compatibility"] = "modern seed stream; not identical to legacy GSL trajectories";
    result["transport_convention"] = "mobility=diffusivity=1 is normalized legacy; independent changes are explicit extensions";
    result["integrator"] = "simultaneous old-state overdamped Euler; additive noise";
    result["observation_phase"] = "current complete state; after all completed steps; initial frame optional";
    result["observation_time"] = "time=absolute_step*dt; stride is relative to sample call; final partial stride included";
    result["array_ownership"] = "independent owned C-contiguous NumPy copies";
    result["axis_order"] = "sample: frame,particle,coordinate; observe: particle,coordinate";
    result["initialization"] = explicit_initial_positions ? "explicitly supplied positions" : "zero positions";
    result["initialization_rng"] = "initialization consumes no RNG";
    result["restart_contract"] = "observations are not restart checkpoints; checkpoint()/from_checkpoint() capture default-component restart state";
    result["topology"] = spec.model == bdmc::Model::Dumbbell
        ? "permanent pairs i<->i+N/2 in initial storage order; bond endpoints are stable IDs"
        : "no fixed bonds";
    if (spec.model == bdmc::Model::Dumbbell) {
        result["spring_coefficient"] = static_cast<double>(spec.dimension) * spec.spring_alpha * spec.spring_alpha;
        result["spring_energy"] = "U=0.5*dimension*spring_alpha^2*sum_pairs(|minimum_image(q)|^2)";
    }
    result["scientific_validation"] = "focused unit/analytic checks only; historical short_test.inp regression deferred";
    return result;
}

py::dict frame_dict(const bdmc::ObservationFrame& frame, const bdmc::SimulationSpec& spec,
                    bool explicit_initial_positions) {
    const std::vector<py::ssize_t> shape{array_size(spec.n_particles), array_size(spec.dimension)};
    py::dict result;
    result["positions"] = copied_array(frame.state.positions, shape);
    result["unwrapped_positions"] = copied_array(frame.state.unwrapped_positions, shape);
    result["forces"] = copied_array(frame.force.values, shape);
    result["potential_energy"] = frame.force.potential_energy;
    result["step"] = frame.state.step;
    result["time"] = frame.state.time;
    result["particle_ids"] = copied_array(frame.state.ids, {array_size(spec.n_particles)});
    result["dynamic_bonds"] = bond_array(frame.dynamic_bonds);
    result["metadata"] = metadata(spec, explicit_initial_positions);
    return result;
}

py::dict batch_dict(const bdmc::ObservationBatch& batch, const bdmc::SimulationSpec& spec,
                    const std::vector<bdmc::ParticleId>& ids, bool explicit_initial_positions) {
    // ML/analysis receives a declared (frame,particle,coordinate) tensor rather
    // than an implicit text-column convention. Ragged dynamic topology stays a
    // per-frame list; padding it with zeros would invent physical bonds.
    const auto frames = array_size(batch.frames.size());
    const auto particles = array_size(spec.n_particles);
    const auto dimension = array_size(spec.dimension);
    const std::vector<py::ssize_t> shape{frames, particles, dimension};
    py::array_t<double> positions(shape), unwrapped(shape), forces(shape);
    py::array_t<double> times(frames), energies(frames);
    py::array_t<std::uint64_t> steps(frames);
    auto* pos = positions.mutable_data();
    auto* unwrap = unwrapped.mutable_data();
    auto* force = forces.mutable_data();
    py::list dynamic;
    for (std::size_t k = 0; k < batch.frames.size(); ++k) {
        const auto& frame = batch.frames[k];
        pos = std::copy(frame.state.positions.begin(), frame.state.positions.end(), pos);
        unwrap = std::copy(frame.state.unwrapped_positions.begin(), frame.state.unwrapped_positions.end(), unwrap);
        force = std::copy(frame.force.values.begin(), frame.force.values.end(), force);
        times.mutable_data()[k] = frame.state.time;
        energies.mutable_data()[k] = frame.force.potential_energy;
        steps.mutable_data()[k] = frame.state.step;
        dynamic.append(bond_array(frame.dynamic_bonds));
    }
    py::dict result;
    result["positions"] = std::move(positions);
    result["unwrapped_positions"] = std::move(unwrapped);
    result["forces"] = std::move(forces);
    result["times"] = std::move(times);
    result["steps"] = std::move(steps);
    result["potential_energy"] = std::move(energies);
    result["particle_ids"] = copied_array(ids, {particles});
    result["fixed_bonds"] = bond_array(batch.fixed_bonds);
    result["dynamic_bonds"] = std::move(dynamic);
    result["metadata"] = metadata(spec, explicit_initial_positions);
    return result;
}

std::vector<double> initial_positions(py::object value, const bdmc::SimulationSpec& spec) {
    if (value.is_none()) return {};
    auto array = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(value);
    if (!array || array.ndim() != 2 || array.shape(0) != array_size(spec.n_particles) ||
        array.shape(1) != array_size(spec.dimension)) {
        throw std::invalid_argument("initial_positions must have shape (n_particles,dimension)");
    }
    return {array.data(), array.data() + array.size()};
}

std::vector<bdmc::ParticleId> initial_ids(py::object value) {
    if (value.is_none()) return {};
    try {
        return value.cast<std::vector<bdmc::ParticleId>>();
    } catch (const py::cast_error&) {
        throw std::invalid_argument("particle_ids must be a sequence of nonnegative 64-bit integers");
    }
}

py::dict neighbor_dict(const bdmc::NeighborResult& result, const bdmc::SpatialSearchSpec& spec) {
    const auto count = array_size(result.pairs.size());
    const auto dimension = array_size(spec.box_lengths.size());
    py::array_t<std::uint64_t> indices({count, py::ssize_t{2}});
    py::array_t<std::uint64_t> ids({count, py::ssize_t{2}});
    py::array_t<double> distances(count);
    py::array_t<double> vectors({count, dimension});
    for (std::size_t k = 0; k < result.pairs.size(); ++k) {
        const auto& pair = result.pairs[k];
        indices.mutable_data()[2*k] = pair.first_index;
        indices.mutable_data()[2*k+1] = pair.second_index;
        ids.mutable_data()[2*k] = pair.first_id;
        ids.mutable_data()[2*k+1] = pair.second_id;
        distances.mutable_data()[k] = pair.distance;
        for (std::size_t c = 0; c < spec.box_lengths.size(); ++c) {
            vectors.mutable_data()[k * spec.box_lengths.size() + c] = pair.displacement[c];
        }
    }
    for (py::array array : std::vector<py::array>{indices, ids, distances, vectors}) {
        array.attr("setflags")(py::arg("write") = false);
    }
    py::dict stats;
    stats["backend"] = bdmc::spatial_backend_name(result.stats.backend);
    stats["n_particles"] = result.stats.n_particles;
    stats["dimension"] = result.stats.dimension;
    stats["occupied_cells"] = result.stats.occupied_cells;
    stats["grid_shape"] = result.stats.grid_shape;
    stats["candidate_pairs"] = result.stats.candidate_pairs;
    stats["pair_count"] = result.pairs.size();
    stats["stored_particle_indices"] = result.stats.stored_particle_indices;
    stats["particle_index_itemsize"] = sizeof(std::size_t);
    stats["cutoff"] = spec.cutoff;
    stats["pair_rule"] = "distinct particles; minimum-image distance <= cutoff; sorted storage indices i<j";
    stats["vector_orientation"] = "minimum_image(position[j]-position[i]); [-L/2,L/2)";
    stats["array_ownership"] = "independent copied C-contiguous arrays; read-only flags";
    stats["storage_contract"] = "sparse occupied cells; streamed candidates; retained pair output O(P); no full-grid allocation";
    if (result.stats.grid_shape.empty()) {
        stats["dense_grid_cells"] = py::none();
    } else {
        // Python arbitrary-precision integer arithmetic reports a theoretical
        // grid even when its product would overflow a native integer type.
        py::object dense_count = py::int_(1);
        for (const auto side : result.stats.grid_shape) dense_count = dense_count.attr("__mul__")(py::int_(side));
        stats["dense_grid_cells"] = dense_count;
    }
    py::dict output;
    output["pair_indices"] = std::move(indices);
    output["pair_ids"] = std::move(ids);
    output["distances"] = std::move(distances);
    output["vectors"] = std::move(vectors);
    output["stats"] = std::move(stats);
    return output;
}

void checkpoint_fields(const py::dict& input, std::initializer_list<const char*> fields,
                        const std::string& context) {
    if (input.size() != fields.size()) {
        throw std::invalid_argument(context + " has missing or unknown fields");
    }
    for (const auto* field : fields) {
        if (!input.contains(py::str(field))) throw std::invalid_argument(context + " is missing " + field);
    }
}

std::uint64_t checkpoint_unsigned(py::handle value, const std::string& field) {
    if (!py::isinstance<py::int_>(value) || py::isinstance<py::bool_>(value)) {
        throw std::invalid_argument(field + " must be a nonnegative integer");
    }
    try {
        return py::cast<std::uint64_t>(value);
    } catch (const py::cast_error&) {
        throw std::invalid_argument(field + " must fit uint64");
    }
}

bool checkpoint_bool(py::handle value, const std::string& field) {
    if (!py::isinstance<py::bool_>(value)) throw std::invalid_argument(field + " must be a boolean");
    return py::cast<bool>(value);
}

void checkpoint_array_dtype(py::handle value, const py::dtype& dtype, const std::string& field) {
    // Ordinary input may be converted to float64, but restart arrays carry an
    // exact schema: reject unintended precision/byte-order changes before copy.
    if (py::isinstance<py::array>(value)) {
        const auto array = py::reinterpret_borrow<py::array>(value);
        if (py::cast<std::string>(array.dtype().attr("str")) != py::cast<std::string>(dtype.attr("str"))) {
            throw std::invalid_argument(field + " has the wrong numeric dtype or byte order");
        }
    }
}

std::vector<double> checkpoint_positions(py::handle value, const bdmc::SimulationSpec& spec,
                                          const std::string& field) {
    checkpoint_array_dtype(value, py::dtype::of<double>(), field);
    auto array = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(value);
    if (!array || array.ndim() != 2 || array.shape(0) != array_size(spec.n_particles) ||
        array.shape(1) != array_size(spec.dimension)) {
        throw std::invalid_argument(field + " must have shape (n_particles,dimension)");
    }
    return {array.data(), array.data() + array.size()};
}

std::vector<bdmc::ParticleId> checkpoint_ids(py::handle value, const bdmc::SimulationSpec& spec) {
    checkpoint_array_dtype(value, py::dtype::of<bdmc::ParticleId>(), "state.particle_ids");
    if (py::isinstance<py::array>(value)) {
        const auto array = py::reinterpret_borrow<py::array>(value);
        if (array.ndim() != 1 || array.shape(0) != array_size(spec.n_particles)) {
            throw std::invalid_argument("state.particle_ids must have shape (n_particles,)");
        }
    }
    auto ids = py::cast<std::vector<bdmc::ParticleId>>(value);
    if (ids.size() != spec.n_particles) throw std::invalid_argument("checkpoint particle ID count mismatch");
    return ids;
}

std::vector<bdmc::Bond> checkpoint_bonds(py::handle value, const std::string& field) {
    checkpoint_array_dtype(value, py::dtype::of<bdmc::ParticleId>(), field);
    if (py::isinstance<py::array>(value)) {
        const auto array = py::reinterpret_borrow<py::array>(value);
        if (array.ndim() != 2 || array.shape(1) != 2) {
            throw std::invalid_argument(field + " must have shape (bonds,2)");
        }
    }
    const auto rows = py::cast<std::vector<std::vector<bdmc::ParticleId>>>(value);
    std::vector<bdmc::Bond> bonds;
    bonds.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.size() != 2) throw std::invalid_argument(field + " must have two endpoints per bond");
        bonds.push_back({row[0], row[1]});
    }
    return bonds;
}

py::dict checkpoint_dict(const bdmc::SimulationCheckpoint& snapshot) {
    py::dict result;
    result["schema"] = snapshot.schema;
    result["version"] = snapshot.version;
    result["core_version"] = snapshot.core_version;
    result["backend"] = snapshot.backend;
    result["spec"] = spec_dict(snapshot.spec);
    py::dict geometry;
    geometry["boundary"] = snapshot.boundary;
    geometry["box_lengths"] = snapshot.box_lengths;
    result["geometry"] = std::move(geometry);
    py::dict state;
    const std::vector<py::ssize_t> shape{array_size(snapshot.spec.n_particles), array_size(snapshot.spec.dimension)};
    state["particle_ids"] = copied_array(snapshot.state.ids, {shape[0]});
    state["positions"] = copied_array(snapshot.state.positions, shape);
    state["unwrapped_positions"] = copied_array(snapshot.state.unwrapped_positions, shape);
    state["step"] = snapshot.state.step;
    state["time"] = snapshot.state.time;
    result["state"] = std::move(state);
    py::dict topology;
    topology["fixed_bonds"] = bond_array(snapshot.fixed_bonds);
    topology["dynamic_bonds"] = bond_array(snapshot.dynamic_bonds);
    result["topology"] = std::move(topology);
    py::dict rng;
    rng["engine"] = snapshot.noise.engine;
    rng["engine_state"] = snapshot.noise.engine_state;
    rng["noise_policy"] = bdmc::noise_name(snapshot.noise.policy);
    rng["has_spare"] = snapshot.noise.has_spare;
    rng["spare"] = snapshot.noise.spare;
    result["rng"] = std::move(rng);
    py::dict initialization;
    initialization["explicit_initial_positions"] = snapshot.explicit_initial_positions;
    result["initialization"] = std::move(initialization);
    return result;
}

bdmc::SimulationCheckpoint parse_checkpoint(const py::dict& input) {
    checkpoint_fields(input, {"schema", "version", "core_version", "backend", "spec", "geometry",
                              "state", "topology", "rng", "initialization"}, "checkpoint");
    bdmc::SimulationCheckpoint snapshot;
    snapshot.schema = py::cast<std::string>(input["schema"]);
    const auto version = checkpoint_unsigned(input["version"], "checkpoint.version");
    if (version != 1) throw std::invalid_argument("unsupported checkpoint schema version");
    snapshot.version = static_cast<std::uint32_t>(version);
    snapshot.core_version = py::cast<std::string>(input["core_version"]);
    snapshot.backend = py::cast<std::string>(input["backend"]);
    const auto spec = py::cast<py::dict>(input["spec"]);
    checkpoint_fields(spec, {"model", "n_particles", "dimension", "box_lengths", "dt", "mobility",
                            "diffusivity", "spring_alpha", "seed", "noise_policy"}, "checkpoint.spec");
    snapshot.spec.model = bdmc::parse_model(py::cast<std::string>(spec["model"]));
    const auto particles = checkpoint_unsigned(spec["n_particles"], "spec.n_particles");
    const auto dimension = checkpoint_unsigned(spec["dimension"], "spec.dimension");
    if (dimension < 1 || dimension > 3) throw std::invalid_argument("checkpoint dimension must be 1, 2 or 3");
    if (particles > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("checkpoint particle count overflows size_t");
    snapshot.spec.n_particles = static_cast<std::size_t>(particles);
    snapshot.spec.dimension = static_cast<std::size_t>(dimension);
    snapshot.spec.box_lengths = py::cast<std::vector<double>>(spec["box_lengths"]);
    snapshot.spec.dt = py::cast<double>(spec["dt"]);
    snapshot.spec.mobility = py::cast<double>(spec["mobility"]);
    snapshot.spec.diffusivity = py::cast<double>(spec["diffusivity"]);
    snapshot.spec.spring_alpha = py::cast<double>(spec["spring_alpha"]);
    snapshot.spec.seed = checkpoint_unsigned(spec["seed"], "spec.seed");
    snapshot.spec.noise_policy = bdmc::parse_noise_policy(py::cast<std::string>(spec["noise_policy"]));
    snapshot.spec.validate();

    const auto geometry = py::cast<py::dict>(input["geometry"]);
    checkpoint_fields(geometry, {"boundary", "box_lengths"}, "checkpoint.geometry");
    snapshot.boundary = py::cast<std::string>(geometry["boundary"]);
    snapshot.box_lengths = py::cast<std::vector<double>>(geometry["box_lengths"]);
    const auto state = py::cast<py::dict>(input["state"]);
    checkpoint_fields(state, {"particle_ids", "positions", "unwrapped_positions", "step", "time"}, "checkpoint.state");
    snapshot.state.dimension = snapshot.spec.dimension;
    snapshot.state.ids = checkpoint_ids(state["particle_ids"], snapshot.spec);
    snapshot.state.positions = checkpoint_positions(state["positions"], snapshot.spec, "state.positions");
    snapshot.state.unwrapped_positions = checkpoint_positions(state["unwrapped_positions"], snapshot.spec, "state.unwrapped_positions");
    snapshot.state.step = checkpoint_unsigned(state["step"], "state.step");
    snapshot.state.time = py::cast<double>(state["time"]);
    const auto topology = py::cast<py::dict>(input["topology"]);
    checkpoint_fields(topology, {"fixed_bonds", "dynamic_bonds"}, "checkpoint.topology");
    snapshot.fixed_bonds = checkpoint_bonds(topology["fixed_bonds"], "topology.fixed_bonds");
    snapshot.dynamic_bonds = checkpoint_bonds(topology["dynamic_bonds"], "topology.dynamic_bonds");
    const auto rng = py::cast<py::dict>(input["rng"]);
    checkpoint_fields(rng, {"engine", "engine_state", "noise_policy", "has_spare", "spare"}, "checkpoint.rng");
    snapshot.noise.engine = py::cast<std::string>(rng["engine"]);
    snapshot.noise.engine_state = py::cast<std::string>(rng["engine_state"]);
    snapshot.noise.policy = bdmc::parse_noise_policy(py::cast<std::string>(rng["noise_policy"]));
    snapshot.noise.has_spare = checkpoint_bool(rng["has_spare"], "rng.has_spare");
    snapshot.noise.spare = py::cast<double>(rng["spare"]);
    const auto initialization = py::cast<py::dict>(input["initialization"]);
    checkpoint_fields(initialization, {"explicit_initial_positions"}, "checkpoint.initialization");
    snapshot.explicit_initial_positions = checkpoint_bool(initialization["explicit_initial_positions"], "initialization.explicit_initial_positions");
    return snapshot;
}

}  // namespace

PYBIND11_MODULE(_core, module) {
    module.doc() = "Composable C++17 Brownian/dumbbell CPU reference; explicit physical and observation contracts.";
    module.attr("__version__") = "0.1.0";
    py::enum_<bdmc::Model>(module, "Model")
        .value("BROWNIAN", bdmc::Model::Brownian)
        .value("DUMBBELL", bdmc::Model::Dumbbell);
    py::enum_<bdmc::NoisePolicy>(module, "NoisePolicy")
        .value("LEGACY_UNIFORM", bdmc::NoisePolicy::LegacyUniform)
        .value("GAUSSIAN", bdmc::NoisePolicy::Gaussian);
    py::enum_<bdmc::SpatialBackend>(module, "SpatialBackend")
        .value("CELL_LIST", bdmc::SpatialBackend::CellList)
        .value("BRUTE_FORCE", bdmc::SpatialBackend::BruteForce);

    module.def("neighbor_pairs", [](py::object positions, std::vector<double> lengths,
                                     double cutoff, const std::string& backend, py::object particle_ids) {
        bdmc::SpatialSearchSpec spec{std::move(lengths), cutoff, bdmc::parse_spatial_backend(backend)};
        spec.validate();
        auto array = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(positions);
        if (!array || array.ndim() != 2 || array.shape(1) != array_size(spec.box_lengths.size())) {
            throw std::invalid_argument("positions must have shape (n_particles,dimension) matching box_lengths");
        }
        const auto n = static_cast<std::size_t>(array.shape(0));
        auto ids = initial_ids(particle_ids);
        if (!particle_ids.is_none() && ids.size() != n) throw std::invalid_argument("spatial particle ID count mismatch");
        std::vector<double> input;
        if (array.size() > 0) input.assign(array.data(), array.data() + array.size());
        bdmc::NeighborResult result;
        {
            py::gil_scoped_release release;
            result = bdmc::SpatialSearch(spec)(input, std::move(ids));
        }
        return neighbor_dict(result, spec);
    }, py::arg("positions"), py::arg("box_lengths"), py::arg("cutoff"),
       py::arg("backend") = "cell_list", py::arg("particle_ids") = py::none(),
       "Return independent read-only neighbor arrays; 0<cutoff<=min(box)/2; inclusive distance threshold.");

    py::class_<bdmc::SimulationSpec>(module, "SimulationSpec")
        .def(py::init([](const std::string& model, std::size_t n_particles,
                         std::size_t dimension, std::vector<double> lengths,
                         double dt, double mobility, double diffusivity,
                         double spring_alpha, std::uint64_t seed,
                         const std::string& noise_policy) {
            bdmc::SimulationSpec spec;
            if (dimension < 1 || dimension > 3) {
                throw std::invalid_argument("dimension must be 1, 2 or 3");
            }
            spec.model = bdmc::parse_model(model);
            spec.n_particles = n_particles;
            spec.dimension = dimension;
            spec.box_lengths = lengths.empty() ? std::vector<double>(dimension, 10.0) : std::move(lengths);
            spec.dt = dt;
            spec.mobility = mobility;
            spec.diffusivity = diffusivity;
            spec.spring_alpha = spring_alpha;
            spec.seed = seed;
            spec.noise_policy = bdmc::parse_noise_policy(noise_policy);
            spec.validate();
            return spec;
        }), py::arg("model") = "brownian", py::arg("n_particles") = 64,
            py::arg("dimension") = 3, py::arg("box_lengths") = std::vector<double>{},
            py::arg("dt") = 0.001, py::arg("mobility") = 1.0,
            py::arg("diffusivity") = 1.0, py::arg("spring_alpha") = 1.0,
            py::arg("seed") = 0, py::arg("noise_policy") = "legacy_uniform")
        .def_property("model", [](const bdmc::SimulationSpec& spec) { return bdmc::model_name(spec.model); },
            [](bdmc::SimulationSpec& spec, const std::string& name) { spec.model = bdmc::parse_model(name); })
        .def_readwrite("n_particles", &bdmc::SimulationSpec::n_particles)
        .def_readwrite("dimension", &bdmc::SimulationSpec::dimension)
        .def_readwrite("box_lengths", &bdmc::SimulationSpec::box_lengths)
        .def_readwrite("dt", &bdmc::SimulationSpec::dt)
        .def_readwrite("mobility", &bdmc::SimulationSpec::mobility)
        .def_readwrite("diffusivity", &bdmc::SimulationSpec::diffusivity)
        .def_readwrite("spring_alpha", &bdmc::SimulationSpec::spring_alpha)
        .def_readwrite("seed", &bdmc::SimulationSpec::seed)
        .def_property("noise_policy", [](const bdmc::SimulationSpec& spec) { return bdmc::noise_name(spec.noise_policy); },
            [](bdmc::SimulationSpec& spec, const std::string& name) { spec.noise_policy = bdmc::parse_noise_policy(name); })
        .def("validate", &bdmc::SimulationSpec::validate)
        .def("to_dict", &spec_dict);

    py::class_<bdmc::Simulation>(module, "Simulation")
        .def(py::init([](bdmc::SimulationSpec spec, py::object positions, py::object ids) {
            // Validate before any array-size multiplication/conversion.
            spec.validate();
            auto particle_ids = initial_ids(ids);
            if (!ids.is_none() && particle_ids.size() != spec.n_particles) {
                throw std::invalid_argument("particle ID count differs from specification");
            }
            return std::make_unique<bdmc::Simulation>(spec, initial_positions(positions, spec), std::move(particle_ids));
        }), py::arg("spec"), py::arg("initial_positions") = py::none(), py::arg("particle_ids") = py::none())
        .def_property_readonly("spec", &bdmc::Simulation::spec)
        .def_property_readonly("metadata", [](const bdmc::Simulation& simulation) {
            return metadata(simulation.spec(), simulation.has_explicit_initial_positions());
        })
        .def_property_readonly("fixed_bonds", [](const bdmc::Simulation& simulation) { return bond_array(simulation.fixed_bonds()); })
        .def_property_readonly("dynamic_bonds", [](const bdmc::Simulation& simulation) { return bond_array(simulation.dynamic_bonds()); })
        .def("advance", &bdmc::Simulation::advance, py::arg("steps"), py::call_guard<py::gil_scoped_release>())
        .def("observe", [](const bdmc::Simulation& simulation) {
            // Native work releases the GIL and uses the simulation's own lock.
            // NumPy/dict construction resumes with the GIL held; never touch
            // Python objects from inside the released native-work scope.
            bdmc::ObservationFrame frame;
            {
                py::gil_scoped_release release;
                frame = simulation.observe();
            }
            auto result = frame_dict(frame, simulation.spec(), simulation.has_explicit_initial_positions());
            result["fixed_bonds"] = bond_array(simulation.fixed_bonds());
            return result;
        })
        .def("sample", [](bdmc::Simulation& simulation, std::uint64_t steps,
                           std::uint64_t stride, bool include_initial) {
            bdmc::ObservationBatch batch;
            {
                py::gil_scoped_release release;
                batch = simulation.sample(steps, stride, include_initial);
            }
            return batch_dict(batch, simulation.spec(), simulation.state().ids,
                              simulation.has_explicit_initial_positions());
        }, py::arg("steps"), py::arg("stride") = 1, py::arg("include_initial") = true)
        .def("checkpoint", [](const bdmc::Simulation& simulation) {
            bdmc::SimulationCheckpoint snapshot;
            {
                py::gil_scoped_release release;
                snapshot = simulation.checkpoint();
            }
            return checkpoint_dict(snapshot);
        })
        .def_static("from_checkpoint", [](const py::dict& input) {
            try {
                const auto snapshot = parse_checkpoint(input);
                py::gil_scoped_release release;
                return bdmc::Simulation::from_checkpoint(snapshot);
            } catch (const py::cast_error&) {
                throw std::invalid_argument("checkpoint contains an invalid field type");
            }
        }, py::arg("snapshot"));

#ifdef BDMC_HAS_EXTRA_BINDINGS
    bind_extra(module);
#endif
}
