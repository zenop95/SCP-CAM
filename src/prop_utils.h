#ifndef AIDA_IO_UTILS_H
#define AIDA_IO_UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <dace/dace.h>
#include <cspice/SpiceUsr.h>
#include <nlohmann/json.hpp>

//
// ================================================================
//                   namespace aida::io
// ================================================================
//

namespace props { namespace io {

// ======================
// TIME UTILITIES
// ======================

long now_ms();
void set_prec(std::ofstream& os, int p = 16);

// ======================
// JSON TOOLS
// ======================

nlohmann::json read_json(const std::string& path);
void write_json(const std::string& path, const nlohmann::json& j);

// Accepts either [a,b,c] or [a] style or [a,b,...] arrays.
// Returns a 1D vector in both cases.
std::vector<double> read_vec(const nlohmann::json& j, const char* key);

// r field readers (scaled position):
//  - (N×M×3) array-of-arrays OR (N×3) when M=1
// Returns flattened list of length N*M, each element size 3.
std::vector<std::vector<double>>
read_r_list(const nlohmann::json& j, const char* key, int N, int M);

// P / e2b (3×3 matrices):
//  - For PoC: list of M matrices OR single 3×3 when M=1
std::vector<std::vector<std::vector<double>>>
read_mat3_list(const nlohmann::json& j, const char* key, int M);

//  - For IPC: N×M list of 3×3 OR N×(3×3) when M=1
std::vector<std::vector<std::vector<double>>>
read_mat3_list(const nlohmann::json& j, const char* key, int N, int M);

// ======================
/* MATRIX SERIALIZATION (C++ → JSON) */
// ======================

// Convert AlgebraicMatrix<double> to row-major vector<vector<double>>
std::vector<std::vector<double>>
to_vec2(const DACE::AlgebraicMatrix<double>& M, int R, int C);

// Convert a 3×3 AlgebraicMatrix<double> to nested vector form
std::vector<std::vector<double>>
to_mat3(const DACE::AlgebraicMatrix<double>& A);

// Convert a list of 3×3 AlgebraicMatrix<double> to nested vectors
std::vector<std::vector<std::vector<double>>>
to_mat3_list(const std::vector<DACE::AlgebraicMatrix<double>>& mats);

// ======================
// DATA STRUCTURES
// ======================

struct AidaParams {
    int flag1{}, flag2{}, flag3{}, gravOrd{};
    double mass{}, A_drag{}, Cd{}, A_srp{}, Cr{};
};

struct InputDataGeneral {
    int N{};
    double Lsc{}, et{}, aMax{};
    std::vector<double> t;
    DACE::AlgebraicMatrix<double> xdum;  // 6×N
    DACE::AlgebraicMatrix<double> udum;  // 3×N
    int order{}, flagRtn{}, flagProp{};
};

struct InputDataDebris {
    int N{};
    double Lsc{}, et{};
    std::vector<double> t;
    DACE::AlgebraicVector<double> x0;  // 6 elements
};

struct InputDataTCA {
    int M{};
    DACE::AlgebraicMatrix<double> xdump; // primary, 6×M
    DACE::AlgebraicMatrix<double> xdums; // secondary, 6×M
};

struct InputDataTrustRegion {
    double dt{};
    double et{};
    DACE::AlgebraicVector<double> x0; // 6
};

// ======================
// JSON PARSERS
// ======================

AidaParams           read_aida_params_json(const std::string& path);
InputDataGeneral     read_input_general_json(const std::string& path);
InputDataDebris      read_input_debris_json(const std::string& path);
InputDataTCA         read_input_tca_json(const std::string& path);
InputDataTrustRegion read_input_trustregion_json(const std::string& path);

// ======================
// SPICE / AIDA
// ======================

void load_spice_kernel(const std::string& kernelPath);
void clear_spice_kernels();

struct AidaCoeff {
    double Bfactor{};
    double SRPC{};
};

AidaCoeff compute_aida_coeff(const AidaParams& p);

}} // namespace aida::io

#endif