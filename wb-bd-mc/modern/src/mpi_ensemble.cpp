#include "bdmc/reference.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Independent complete reference runs are the unit of MPI ownership. A rank
// never writes a file, shares a random stream, or decomposes a physical system.
namespace {
constexpr std::size_t max_manifest_bytes = 1024 * 1024;
constexpr std::size_t max_output_bytes = 4 * 1024 * 1024;
constexpr std::size_t max_jobs = 256;
const char* columns = "run_id\tcondition_id\tlineage_id\tmodel\tseed\tparameter\tcount\tsteps\tdt\tdiffusion\tmobility\tbox_length\tsigma\tq0\tshear_rate\twarmup";

void require(bool condition, const std::string& message) {
    if (!condition) throw std::invalid_argument(message);
}
std::string quoted(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32 || c > 126) out << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << unsigned(c) << std::dec;
        else out << c;
    }
    out << '"';
    return out.str();
}
std::uint64_t integer(const std::string& value, const std::string& name, std::uint64_t low,
                      std::uint64_t high) {
    std::uint64_t result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    require(parsed.ec == std::errc{} && parsed.ptr == value.data()+value.size() && result >= low && result <= high,
            name + " is outside its integer bounds");
    return result;
}
double number(const std::string& value, const std::string& name, double low, double high,
              bool positive = false) {
    std::istringstream stream(value);
    stream.imbue(std::locale::classic());
    double result = 0;
    stream >> std::noskipws >> result;
    require(stream && stream.eof() && std::isfinite(result) && result >= low && result <= high && (!positive || result > 0),
            name + " is outside its finite numeric bounds");
    return result;
}
void identity(const std::string& value) {
    require(!value.empty() && value.size() <= 128, "identities must have 1..128 ASCII characters");
    for (unsigned char c : value)
        require((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ||
                c == '_' || c == '-' || c == '.' || c == ':', "identity has an unsupported character");
}
std::vector<std::string> split(const std::string& value) {
    std::vector<std::string> result;
    std::size_t start = 0;
    while (true) {
        const auto end = value.find('\t', start);
        result.push_back(value.substr(start, end == std::string::npos ? end : end-start));
        if (end == std::string::npos) break;
        start = end+1;
    }
    return result;
}
struct Job {
    std::string run, condition, lineage, model;
    std::uint64_t seed, steps, warmup;
    std::size_t count;
    double parameter, dt, diffusion, mobility, box, sigma, q0, rate;
};
std::vector<Job> parse_jobs(const std::string& manifest) {
    require(!manifest.empty() && manifest.size() <= max_manifest_bytes, "manifest must be 1..1048576 bytes");
    std::istringstream input(manifest);
    std::string line;
    require(bool(std::getline(input,line)) && line == "bdmc.mpi-jobs.v1", "unsupported MPI job manifest schema");
    require(bool(std::getline(input,line)) && line == columns, "MPI manifest column header does not match schema");
    std::vector<Job> jobs;
    std::set<std::string> seen;
    std::set<std::uint64_t> seeds;
    long double work = 0;
    while (std::getline(input,line)) {
        require(jobs.size() < max_jobs, "MPI manifest has more than 256 jobs");
        const auto field = split(line);
        require(field.size() == 16, "each MPI manifest row requires exactly 16 tab-separated fields");
        Job job{};
        job.run=field[0]; job.condition=field[1]; job.lineage=field[2]; job.model=field[3];
        identity(job.run); identity(job.condition); identity(job.lineage);
        require(seen.insert(job.run).second,"duplicate run_id in MPI manifest");
        require(job.model == "repulsive" || job.model == "hookean" || job.model == "fene", "unsupported model in MPI manifest");
        job.seed=integer(field[4],"seed",0,std::numeric_limits<std::uint64_t>::max());
        require(seeds.insert(job.seed).second,"duplicate simulation seed in independent-run MPI manifest");
        job.parameter=number(field[5],"parameter",0,1e6,job.model != "repulsive");
        job.count=integer(field[6],"count",1,job.model == "repulsive" ? 512 : 4096);
        job.steps=integer(field[7],"steps",1,100000);
        job.dt=number(field[8],"dt",0,1,true);
        job.diffusion=number(field[9],"diffusion",0,1e6);
        job.mobility=number(field[10],"mobility",0,1e6,true);
        job.box=number(field[11],"box_length",0,1e6,true);
        job.sigma=number(field[12],"sigma",0,1e6,true);
        job.q0=number(field[13],"q0",0,1e6,true);
        job.rate=number(field[14],"shear_rate",-1e6,1e6);
        job.warmup=integer(field[15],"warmup",0,100000);
        require(job.model != "repulsive" || (job.sigma <= job.box/2 && job.rate == 0),
                "repulsive jobs require sigma<=box_length/2 and no shear");
        require(job.model == "repulsive" || job.mobility*job.parameter*job.dt < 1,
                "dumbbell jobs require mobility*spring_constant*dt<1");
        const long double count = job.count;
        work += (job.steps+job.warmup)*(job.model == "repulsive" ? count*count : 2*count);
        require(work <= 300000000.L,"MPI batch exceeds its 300 million update-work bound");
        jobs.push_back(job);
    }
    require(!jobs.empty(),"MPI manifest contains no jobs");
    return jobs;
}

struct Options { std::string path, backend="serial"; int threads=1; };
Options options(int argc, char** argv) {
    Options result;
    std::set<std::string> seen;
    for (int i=1; i<argc; i+=2) {
        const std::string key=argv[i];
        require(i+1<argc && seen.insert(key).second,"options require distinct name/value pairs");
        const std::string value=argv[i+1];
        if (key == "--jobs-file") result.path=value;
        else if (key == "--backend") result.backend=value;
        else if (key == "--threads") result.threads=static_cast<int>(integer(value,"threads",1,256));
        else throw std::invalid_argument("unknown option: "+key);
    }
    require(!result.path.empty(),"usage: bdmc_mpi_ensemble --jobs-file PATH [--backend serial|openmp|cuda|metal|accelerate] [--threads N]");
    require(result.backend == "serial" || result.backend == "openmp" || result.backend == "cuda" ||
            result.backend == "metal" || result.backend == "accelerate", "unsupported execution backend");
    return result;
}
std::string read_manifest(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    require(input.good(),"cannot open MPI job manifest");
    std::string value(max_manifest_bytes+1,'\0');
    input.read(&value[0],static_cast<std::streamsize>(value.size()));
    value.resize(static_cast<std::size_t>(input.gcount()));
    require(!input.bad(),"failed reading MPI job manifest");
    (void)parse_jobs(value);
    return value;
}
void broadcast(std::string& value, int rank) {
    int length = rank == 0 ? static_cast<int>(value.size()) : 0;
    MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
    if (rank != 0) value.resize(static_cast<std::size_t>(length));
    MPI_Bcast(value.data(),length,MPI_CHAR,0,MPI_COMM_WORLD);
}

// Canonical little-endian binary64 FNV-1a is a comparison fingerprint, not a
// cryptographic identity. Python records SHA256 of the executable and manifest.
std::string fingerprint(const std::vector<double>& values) {
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                  "MPI fingerprints require IEEE binary64");
    std::uint64_t hash=14695981039346656037ULL;
    for (double value : values) {
        std::uint64_t bits;
        std::memcpy(&bits,&value,sizeof(bits));
        for (unsigned int byte=0; byte<8; ++byte) { hash ^= (bits>>(8*byte))&255; hash *= 1099511628211ULL; }
    }
    std::ostringstream out; out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}
template<class Values> void array(std::ostream& out, const Values& values) {
    out << '[';
    for (std::size_t i=0; i<values.size(); ++i) { if(i) out << ','; out << values[i]; }
    out << ']';
}
std::string execute(const Job& job, const Options& option) {
    std::ostringstream out;
    out.imbue(std::locale::classic()); out << std::setprecision(17);
    out << "{\"run_id\":" << quoted(job.run) << ",\"condition_id\":" << quoted(job.condition)
        << ",\"lineage_id\":" << quoted(job.lineage) << ",\"model\":" << quoted(job.model)
        << ",\"seed\":" << job.seed << ",\"status\":\"success\",\"force_search\":"
        << quoted(job.model == "repulsive" ? (option.backend == "serial" ? "cell_list" : "brute_force") : "not_applicable")
        << ",\"spec\":{\"parameter\":" << job.parameter
        << ",\"count\":" << job.count << ",\"steps\":" << job.steps << ",\"dt\":" << job.dt
        << ",\"diffusion\":" << job.diffusion << ",\"mobility\":" << job.mobility
        << ",\"box_length\":" << job.box << ",\"sigma\":" << job.sigma << ",\"q0\":" << job.q0
        << ",\"shear_rate\":" << job.rate << ",\"warmup\":" << job.warmup << "},\"summary\":{";
    if (job.model == "repulsive") {
        bdmc::RepulsiveSpec spec;
        spec.n_particles=job.count; spec.box_length=job.box; spec.dt=job.dt; spec.diffusion=job.diffusion;
        spec.mobility=job.mobility; spec.epsilon=job.parameter; spec.sigma=job.sigma; spec.seed=job.seed;
        spec.execution_backend=option.backend; spec.backend_threads=option.threads;
        // Nonserial force kernels explicitly implement all pairs. Do not pass
        // the serial default cell-list contract into their validated specs.
        if(option.backend != "serial") spec.force_search=bdmc::SpatialBackend::BruteForce;
        bdmc::RepulsiveSimulation sim(spec);
        const auto initial=sim.observe();
        sim.advance(job.warmup); const auto origin=sim.observe();
        sim.advance(job.steps); const auto frame=sim.observe();
        double msd=0;
        for (std::size_t k=0; k<frame.unwrapped_positions.size(); ++k) {
            const double delta=frame.unwrapped_positions[k]-origin.unwrapped_positions[k];
            msd+=delta*delta/static_cast<double>(job.count);
        }
        require(std::isfinite(msd),"MSD overflow");
        out << "\"potential_energy_per_particle\":" << frame.potential_energy/job.count
            << ",\"finite_lag_msd\":" << msd << ",\"pair_virial\":"; array(out,frame.pair_virial);
        out << "},\"initial_state_fnv1a64\":" << quoted(fingerprint(initial.unwrapped_positions))
            << ",\"origin_state_fnv1a64\":" << quoted(fingerprint(origin.unwrapped_positions))
            << ",\"final_state_fnv1a64\":" << quoted(fingerprint(frame.unwrapped_positions));
    } else {
        bdmc::DumbbellSpec spec;
        spec.n_dumbbells=job.count; spec.dt=job.dt; spec.diffusion=job.diffusion; spec.mobility=job.mobility;
        spec.spring_constant=job.parameter; spec.spring_law=bdmc::parse_spring_law(job.model);
        spec.q0=job.q0; spec.protocol=bdmc::FlowProtocol::SimpleShear; spec.shear_rate=job.rate; spec.seed=job.seed;
        spec.execution_backend=option.backend; spec.backend_threads=option.threads;
        bdmc::DumbbellSimulation sim(spec);
        const auto initial=sim.observe();
        sim.advance(job.warmup); const auto origin=sim.observe();
        sim.advance(job.steps); const auto frame=sim.observe();
        out << "\"connector_second_moment\":" << frame.conformation[0]+frame.conformation[4]+frame.conformation[8]
            << ",\"polymer_shear_stress_per_dumbbell\":" << frame.polymer_stress[1]
            << ",\"first_normal_stress_per_dumbbell\":" << frame.polymer_stress[0]-frame.polymer_stress[4]
            << ",\"mean_energy\":" << frame.mean_energy << ",\"conformation\":"; array(out,frame.conformation);
        out << ",\"polymer_stress\":"; array(out,frame.polymer_stress);
        out << "},\"initial_state_fnv1a64\":" << quoted(fingerprint(initial.connectors))
            << ",\"origin_state_fnv1a64\":" << quoted(fingerprint(origin.connectors))
            << ",\"final_state_fnv1a64\":" << quoted(fingerprint(frame.connectors));
    }
    out << ",\"step\":" << job.steps+job.warmup << ",\"time\":" << (job.steps+job.warmup)*job.dt
        << ",\"observation_lag\":" << job.steps*job.dt << '}';
    return out.str();
}
} // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    Options option;
    std::string manifest, setup_error;
    if(rank == 0) {
        try { option=options(argc,argv); manifest=read_manifest(option.path); }
        catch (const std::exception& error) { setup_error=error.what(); }
    }
    broadcast(setup_error,rank);
    if(!setup_error.empty()) {
        if(rank == 0) std::cout << "{\"schema\":\"bdmc.mpi-ensemble.v1\",\"status\":\"failed\",\"phase\":\"input\",\"error\":" << quoted(setup_error) << "}\n";
        MPI_Finalize(); return 2;
    }
    broadcast(manifest,rank); broadcast(option.backend,rank);
    MPI_Bcast(&option.threads,1,MPI_INT,0,MPI_COMM_WORLD);
    const auto jobs=parse_jobs(manifest); // Identical validated bounded bytes on every rank.
    std::string local;
    int local_failure=0;
    for (std::size_t index=static_cast<std::size_t>(rank); index<jobs.size(); index+=static_cast<std::size_t>(size)) {
        try { local+=std::to_string(index)+'\t'+execute(jobs[index],option)+'\n'; }
        catch (const std::exception& error) {
            local_failure=1;
            // Keep the original seed and identity. Other independent jobs may
            // finish, but the overall batch fails and the adapter raises.
            local+=std::to_string(index)+"\t{\"run_id\":"+quoted(jobs[index].run)+",\"condition_id\":"+quoted(jobs[index].condition)
                +",\"lineage_id\":"+quoted(jobs[index].lineage)+",\"model\":"+quoted(jobs[index].model)
                +",\"seed\":"+std::to_string(jobs[index].seed)+",\"status\":\"failed\",\"error\":"+quoted(std::string(error.what()).substr(0,2048))+"}\n";
        }
    }
    int failed=0;
    MPI_Allreduce(&local_failure,&failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    int length=static_cast<int>(local.size());
    std::vector<int> lengths(rank == 0 ? static_cast<std::size_t>(size) : 0);
    MPI_Gather(&length,1,MPI_INT,lengths.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    std::vector<int> offsets;
    std::vector<char> received;
    int total=0, gather_error=0;
    if(rank == 0) {
        offsets.resize(lengths.size());
        for(std::size_t i=0;i<lengths.size();++i) { offsets[i]=total; total+=lengths[i]; }
        if(total<0 || static_cast<std::size_t>(total)>max_output_bytes) gather_error=1;
        else received.resize(static_cast<std::size_t>(total));
    }
    MPI_Bcast(&gather_error,1,MPI_INT,0,MPI_COMM_WORLD);
    if(gather_error) {
        if(rank == 0) std::cout << "{\"schema\":\"bdmc.mpi-ensemble.v1\",\"status\":\"failed\",\"error\":\"MPI result exceeds 4MiB bound\"}\n";
        MPI_Finalize(); return 3;
    }
    MPI_Gatherv(local.data(),length,MPI_CHAR,received.data(),lengths.data(),offsets.data(),MPI_CHAR,0,MPI_COMM_WORLD);
    if(rank == 0) {
        std::vector<std::string> ordered(jobs.size());
        std::istringstream collected(std::string(received.begin(),received.end()));
        std::string line;
        while(std::getline(collected,line)) {
            const auto tab=line.find('\t');
            const auto index=integer(line.substr(0,tab),"result index",0,jobs.size()-1);
            ordered[index]=line.substr(tab+1);
        }
        std::cout << "{\"schema\":\"bdmc.mpi-ensemble.v1\",\"status\":" << quoted(failed ? "failed" : "success")
            << ",\"execution\":{\"distribution\":\"independent_runs_round_robin\",\"ranks\":" << size
            << ",\"backend\":" << quoted(option.backend) << ",\"threads_per_rank\":" << option.threads
            << "},\"noise_policy\":\"gaussian\",\"failure_policy\":\"no retries, no replacement seeds, unsuccessful batch fails\",\"results\":[";
        for(std::size_t i=0;i<ordered.size();++i) { if(i) std::cout << ','; std::cout << ordered[i]; }
        std::cout << "]}\n";
    }
    MPI_Finalize(); return failed ? 3 : 0;
}
