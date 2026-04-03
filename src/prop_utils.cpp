#include "prop_utils.h"

namespace props { namespace io {

using namespace std;
using namespace DACE;
using namespace std::chrono;

// =====================================================
// TIME UTILITIES
// =====================================================

long now_ms() {
    return time_point_cast<milliseconds>(system_clock::now())
        .time_since_epoch().count();
}

void set_prec(std::ofstream& os, int p) {
    os << std::setprecision(p);
}

// =====================================================
// JSON UTILITIES
// =====================================================

nlohmann::json read_json(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open JSON file: " + path);
    nlohmann::json j; f >> j;
    return j;
}

void write_json(const std::string& path, const nlohmann::json& j) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write JSON file: " + path);
    f << j.dump(2);
}


std::vector<double> read_vec(const nlohmann::json& j, const char* key) {
    const auto& x = j.at(key);
    std::vector<double> out;
    if (x.is_array()) {
        out.reserve(x.size());
        for (const auto& e : x) out.push_back(e.get<double>());
    } else if (x.is_number()) {
        out.push_back(x.get<double>());
    } else {
        throw std::runtime_error(std::string("Field '") + key + "' must be number or array");
    }
    return out;
}

std::vector<std::vector<double>>
read_r_list(const nlohmann::json& j, const char* key, int N, int M) {
    const auto& x = j.at(key);
    if (!x.is_array())
        throw std::runtime_error(std::string("Field '") + key + "' must be array");

    std::vector<std::vector<double>> out;
    out.reserve(N * M);

    // Case: M==1 and x is (N×3)
    if (M == 1 &&
        x.size() == static_cast<size_t>(N) &&
        x[0].is_array() && x[0].size() == 3 && x[0][0].is_number()) {
        for (int i=0;i<N;i++)
            out.push_back(x[i].get<std::vector<double>>());
        return out;
    }

    // Case: (N×M×3)
    if (x.size() == static_cast<size_t>(N) &&
        x[0].is_array() && !x[0].empty() &&
        x[0][0].is_array() && x[0][0].size() == 3) {
        for (int i=0;i<N;i++)
            for (int m=0;m<M;m++)
                out.push_back(x[i][m].get<std::vector<double>>());
        return out;
    }

    throw std::runtime_error(std::string("Field '") + key + "' has unexpected shape for r");
}

std::vector<std::vector<std::vector<double>>>
read_mat3_list(const nlohmann::json& j, const char* key, int M) {
    const auto& x = j.at(key);
    if (!x.is_array())
        throw std::runtime_error(std::string("Field '") + key + "' must be array");

    std::vector<std::vector<std::vector<double>>> out;

    // Case A: list of 3x3 matrices (size M)
    if (!x.empty() && x[0].is_array() && x[0].size() == 3 && x[0][0].is_array()) {
        out.reserve(x.size());
        for (const auto& m : x) out.push_back(m.get<std::vector<std::vector<double>>>());
        if (static_cast<int>(out.size()) != M && M != 1) {
            throw std::runtime_error(std::string("'") + key + "' list count != M");
        }
        return out;
    }

    // Case B: single 3x3 matrix when M==1
    if (M == 1 &&
        x.size() == 3 && x[0].is_array() && x[0].size() == 3 &&
        x[0][0].is_number()) {
        out.push_back(x.get<std::vector<std::vector<double>>>());
        return out;
    }

    throw std::runtime_error(std::string("Field '") + key + "' has unexpected structure (M variant)");
}

std::vector<std::vector<std::vector<double>>>
read_mat3_list(const nlohmann::json& j, const char* key, int N, int M) {
    const auto& x = j.at(key);
    if (!x.is_array())
        throw std::runtime_error(std::string("Field '") + key + "' must be array");

    std::vector<std::vector<std::vector<double>>> out;
    out.reserve(N * M);

    // Case A: M==1 and x is (N×3×3)
    if (M == 1 &&
        x.size() == static_cast<size_t>(N) &&
        x[0].is_array() && x[0].size() == 3 && x[0][0].is_array() &&
        x[0][0].size() == 3 && x[0][0][0].is_number()) {
        for (int i=0;i<N;i++)
            out.push_back(x[i].get<std::vector<std::vector<double>>>());
        return out;
    }

    // Case B: (N×M×3×3)
    if (x.size() == static_cast<size_t>(N) &&
        x[0].is_array() && !x[0].empty() &&
        x[0][0].is_array() && x[0][0].size() == 3 &&
        x[0][0][0].is_array() && x[0][0][0].size() == 3) {
        for (int i=0;i<N;i++)
            for (int m=0;m<M;m++)
                out.push_back(x[i][m].get<std::vector<std::vector<double>>>());
        return out;
    }

    throw std::runtime_error(std::string("Field '") + key + "' has unexpected structure (N,M variant)");
}

// =====================================================
// AIDA PARAMS (JSON)
// =====================================================

AidaParams read_aida_params_json(const std::string& path) {
    auto j = read_json(path);

    AidaParams p;
    p.flag1   = j.at("flag1").get<int>();
    p.flag2   = j.at("flag2").get<int>();
    p.flag3   = j.at("flag3").get<int>();
    p.gravOrd = j.at("gravOrd").get<int>();

    p.mass    = j.at("mass").get<double>();
    p.A_drag  = j.at("A_drag").get<double>();
    p.Cd      = j.at("Cd").get<double>();
    p.A_srp   = j.at("A_srp").get<double>();
    p.Cr      = j.at("Cr").get<double>();

    return p;
}

// =====================================================
// propAida JSON INPUT
// =====================================================

InputDataGeneral read_input_general_json(const std::string& path)
{
    auto j = read_json(path);

    InputDataGeneral d;

    d.N       = j.at("N").get<int>();
    d.t       = j.at("t").get<std::vector<double>>();
    d.et      = j.at("et").get<double>();
    d.order   = j.at("order").get<int>();
    d.aMax    = j.at("uMax").get<double>();
    d.Lsc     = j.at("scaling").get<double>();
    d.flagRtn = j.at("flagRtn").get<int>();
    d.flagProp= j.at("flagProp").get<int>();

    // udum
    auto u_json = j.at("u").get<std::vector<std::vector<double>>>();
    d.udum = AlgebraicMatrix<double>(3, d.N);
    for (int i = 0; i < d.N; i++)
        for (int k = 0; k < 3; k++)
            d.udum.at(k,i) = u_json[i][k];

    // xdum
    d.xdum = AlgebraicMatrix<double>(6, d.N);

    if (d.flagProp == 1) {
        auto x0 = j.at("x0").get<std::vector<double>>();
        for (int k = 0; k < 6; k++)
            d.xdum.at(k,0) = x0[k];
    }
    else {
        auto traj = j.at("oldTraj").get<std::vector<std::vector<double>>>();
        for (int i = 0; i < d.N; i++)
            for (int k = 0; k < 6; k++)
                d.xdum.at(k,i) = traj[i][k];
    }

    return d;
}

// =====================================================
// propAidaDebris JSON INPUT
// =====================================================

InputDataDebris read_input_debris_json(const std::string& path)
{
    auto j = read_json(path);

    InputDataDebris d;
    d.N   = j.at("N").get<int>();
    d.Lsc = j.at("scaling").get<double>();
    d.et  = j.at("et").get<double>();
    d.t   = j.at("t").get<std::vector<double>>();

    auto x0v = j.at("x0").get<std::vector<double>>();
    d.x0 = AlgebraicVector<double>(6);
    for (int i = 0; i < 6; i++)
        d.x0[i] = x0v[i];

    return d;
}

// =====================================================
// TCA JSON INPUT
// =====================================================

InputDataTCA read_input_tca_json(const std::string& path)
{
    auto j = read_json(path);

    InputDataTCA d;
    d.M = j.at("M").get<int>();

    auto prim  = j.at("xdump").get<std::vector<std::vector<double>>>();
    auto sec   = j.at("xdums").get<std::vector<std::vector<double>>>();

    d.xdump = AlgebraicMatrix<double>(6, d.M);
    d.xdums = AlgebraicMatrix<double>(6, d.M);

    for (int i = 0; i < d.M; ++i) {
        for (int k = 0; k < 6; ++k)
            d.xdump.at(k,i) = prim[i][k];
        for (int k = 0; k < 6; ++k)
            d.xdums.at(k,i) = sec[i][k];
    }

    return d;
}

// =====================================================
// Trust Region JSON INPUT
// =====================================================

InputDataTrustRegion read_input_trustregion_json(const std::string& path)
{
    auto j = read_json(path);

    InputDataTrustRegion d;
    d.dt = j.at("dt").get<double>();
    auto x0v = j.at("x0").get<std::vector<double>>();
    d.x0 = AlgebraicVector<double>(6);
    for (int i = 0; i < 6; ++i)
        d.x0[i] = x0v[i];

    return d;
}

// =====================================================
// SPICE
// =====================================================

void load_spice_kernel(const std::string& kernelPath) {
    furnsh_c(kernelPath.c_str());
}

void clear_spice_kernels() {
    kclear_c();
}

// =====================================================
// Compute Bfactor / SRPC
// =====================================================

AidaCoeff compute_aida_coeff(const AidaParams& p)
{
    AidaCoeff c;
    c.Bfactor = p.Cd * p.A_drag / p.mass;
    c.SRPC    = p.Cr * p.A_srp  / p.mass;
    return c;
}

}}