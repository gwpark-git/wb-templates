#include "file_IO.h"
#include <cmath>
#include <set>
#include <filesystem>
namespace {
std::string trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}
}
COND::COND(const char* filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error(std::string("cannot open configuration: ") + filename);
  std::string line;
  long number = 0;
  while (std::getline(in, line)) {
    ++number;
    line = trim(line.substr(0, line.find('#')));
    if (line.empty()) continue;
    const auto equal = line.find('=');
    if (equal == std::string::npos) throw std::runtime_error("expected key=value at line " + std::to_string(number));
    auto key = trim(line.substr(0, equal)), value = trim(line.substr(equal+1));
    if (key.empty() || value.empty() || !values_.emplace(key, value).second)
      throw std::runtime_error("empty or duplicate configuration field at line " + std::to_string(number));
  }
  // These controls were effectively optional in old input files. Defaults are
  // explicit here; numerical and physical inputs still have to be supplied.
  values_.emplace("inlet_BC", "pressure");
  values_.emplace("CONT", "FALSE");
  values_.emplace("DETAILED_REPORT", "FALSE");
  values_.emplace("ACCOUNT_CAKE_THICKNESS", "FALSE");
  validate();
}
const std::string& COND::operator()(const std::string& key) const {
  auto it = values_.find(key);
  if (it == values_.end()) throw std::runtime_error("missing configuration field: " + key);
  return it->second;
}
void COND::validate() const {
  const std::set<std::string> numeric = {"N_iter","weight","Nz","Nr","R_channel","L_channel","Lp","eta0","T","ref_Pout","ref_DTP","phi_bulk","a_particle","phi_freeze","phiw_init","TOL_phiw","TOL_delta_f","gamma"};
  auto value = [&](const std::string& key) {
    const auto& raw = (*this)(key);
    std::size_t parsed = 0;
    double result = std::stod(raw, &parsed);
    if (parsed != raw.size() || !std::isfinite(result)) throw std::runtime_error("invalid finite number: " + key);
    return result;
  };
  for (const auto& key : numeric) {
    const double v = value(key);
    if (key != "ref_Pout" && v <= 0) throw std::runtime_error("expected positive value: " + key);
  }
  for (const auto& key : {"N_iter", "Nz", "Nr"}) {
    if (value(key) != std::floor(value(key))) throw std::runtime_error("expected integer: " + std::string(key));
  }
  if (value("Nz") < 3 || value("Nr") < 3) throw std::runtime_error("Nz and Nr must be at least 3");
  // A destination can raise these bounds deliberately after its resource review.
  if (value("Nz") > 400 || value("Nr") > 4000 || value("N_iter") > 100000)
    throw std::runtime_error("configuration exceeds template resource bounds");
  if (value("weight") > 1 || value("gamma") > 1 || value("phi_freeze") >= 1 ||
      value("phi_bulk") >= value("phi_freeze") || value("phiw_init") >= 1)
    throw std::runtime_error("invalid relaxation, hydrodynamic radius ratio or concentration domain");
  const auto& bc = (*this)("inlet_BC");
  if (bc != "pressure" && bc != "velocity") throw std::runtime_error("inlet_BC must be pressure or velocity");
  if (value(bc == "pressure" ? "DLP" : "u_inlet") <= 0) throw std::runtime_error("inlet scale must be positive");
  for (const auto& flag : {"CONT", "DETAILED_REPORT", "ACCOUNT_CAKE_THICKNESS"}) {
    if ((*this)(flag) != "TRUE" && (*this)(flag) != "FALSE") throw std::runtime_error("expected TRUE/FALSE: " + std::string(flag));
  }
  if ((*this)("CONT") == "TRUE") {
    if (!std::filesystem::is_regular_file((*this)("FN_CONT"))) throw std::runtime_error("continuation file unavailable");
  }
  if ((*this)("ACCOUNT_CAKE_THICKNESS") == "TRUE" && value("membrane_thickness") <= 0)
    throw std::runtime_error("membrane_thickness must be positive");
  const std::map<std::string, std::set<std::string>> laws = {
    {"fcn_Dc_div_D0", {"SPHS","HS","SID","CT"}},
    {"fcn_eta_div_eta0", {"SPHS","HS","TEST","TEST_EQ","CT"}},
    {"fcn_Pi", {"CS","NO"}}};
  for (const auto& law : laws) {
    const auto& selected = (*this)(law.first);
    if (!law.second.count(selected) && !std::filesystem::is_regular_file(selected))
      throw std::runtime_error("unknown law or unavailable table: " + selected);
  }
  (*this)("fout");
  const std::set<std::string> controls = {"inlet_BC","DLP","u_inlet","CONT","FN_CONT","DETAILED_REPORT","ACCOUNT_CAKE_THICKNESS","membrane_thickness","fcn_Dc_div_D0","fcn_eta_div_eta0","fcn_Pi","fout","phi_cake","porosity"};
  for (const auto& kv : values_) {
    if (!numeric.count(kv.first) && !controls.count(kv.first)) throw std::runtime_error("unknown configuration field: " + kv.first);
  }
  // The original PARAMETER_SET hardcodes phi_cake=0.64 and never uses porosity.
  // A warning preserves that behavior without silently adopting an input value.
  if (values_.count("phi_cake") || values_.count("porosity"))
    std::cerr << "WARNING: legacy phi_cake/porosity fields are ignored; phi_cake remains 0.64\n";
}
int COND::cond_print() const {
  for (const auto& kv : values_) std::cout << kv.first << '=' << kv.second << '\n';
  return 0;
}

void finish_output(std::ofstream& output, const std::string& path) {
  output.flush();
  if (!output) throw std::runtime_error("output write failed: " + path);
  output.close();
  if (!output) throw std::runtime_error("output close failed: " + path);
}
