#include <dace/dace.h>
#include <nlohmann/json.hpp>
#include <chrono>
#include "prop_utils.h"
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"
#include "dynorb/AIDA.h"
#include "dynorb/AIDAwrappers.h"

using json = nlohmann::json;
using namespace std;
using namespace DACE;
using namespace props::io;

std::vector<std::vector<double>>
to_vec2m(const AlgebraicMatrix<double>& M, int R, int C) {
    std::vector<std::vector<double>> out(R, std::vector<double>(C));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
            out[r][c] = M.at(r, c);
    return out;
}

int main() try {

    long t0 = now_ms();

    // ------------------------------------------------------------
    // READ INPUT JSON
    // ------------------------------------------------------------
    json jin = read_json("./data_sharing/initial_state.json");

    int N     = jin["N"].get<int>();
    double Lsc = jin["scaling"].get<double>();
    double dt  = jin["dt"].get<double>();
    double et0 = jin["et0"].get<double>();
    int flagProp = jin["flagProp"].get<int>();     // 1=init prop, 2=use newTraj

    // Read initial state trajectory X0
    AlgebraicMatrix<double> X0(6, N);
    const auto& x0_json = jin["x0"];
    if (!x0_json.is_array())
        throw runtime_error("x0 must be array");

    if (flagProp == 1) {
        // x0_json is a length‑6 vector
        for (int j=0;j<6;j++) X0.at(j,0) = x0_json[j].get<double>();
    } else {
        // x0_json is N x 6
        if ((int)x0_json.size() != N)
            throw runtime_error("x0 array length != N");
        for (int i=0;i<N;i++)
            for (int j=0;j<6;j++)
                X0.at(j,i) = x0_json[i][j].get<double>();
    }

    // ------------------------------------------------------------
    // READ AIDA INIT JSON
    // ------------------------------------------------------------
    AidaParams ap = read_aida_params_json("./data_sharing/AIDA_init.json");
    AidaCoeff coeff = compute_aida_coeff(ap);
    int flags[3] = { ap.flag1, ap.flag2, ap.flag3 };

    // ------------------------------------------------------------
    // DA INIT
    // ------------------------------------------------------------
    DA::init(2,6);
    DA::setEps(1e-15);

    load_spice_kernel("./data/kernel.txt");

    AIDAScaledDynamics<DA> dyn(
        "./data/gravmodels/egm2008",
        ap.gravOrd,
        flags,
        coeff.Bfactor,
        coeff.SRPC
    );

    // ------------------------------------------------------------
    // OUTPUT BUFFERS
    // ------------------------------------------------------------
    vector<array<double,6>> state;        // N × 6
    vector<array<double,2>> llConst;      // N × 2
    vector<vector<vector<double>>> dynMaps; // N × (6x6)
    vector<vector<vector<double>>> llMaps;  // N × (2x3)
    vector<array<double,6>> nuOut;          // N × 6

    // DA vectors
    AlgebraicVector<DA> x0(6), x(6), u(3);
    u = {0.0,0.0,0.0};

    // Initialize x0 for node 0:
    for (int j=0;j<6;j++)
        x0[j] = X0.at(j,0) + DA(j+1);

    // ------------------------------------------------------------
    // NODE 0 OUTPUT
    // ------------------------------------------------------------
    {
        array<double,6> s{};
        for (int j=0;j<6;j++) s[j] = cons(x0[j]);
        state.push_back(s);

        auto STM = astro::stmDace(x0, 6,6);
        dynMaps.push_back(to_vec2m(STM,6,6));

        auto LLA = astro::llaMaps(x0, et0, Lsc);
        auto LLAmap = astro::stmDace(LLA, 6,3); // 6x3
        // BUT your MATLAB expects 2×3 maps → use only LAT,LON rows 0,1
        vector<vector<double>> ll(2, vector<double>(3));
        for (int r=0;r<2;r++)
            for (int c=0;c<3;c++)
                ll[r][c] = LLAmap.at(r,c);
        llMaps.push_back(ll);

        llConst.push_back({ cons(LLA[0]), cons(LLA[1]) });  // 2×1 (lat, lon)

        auto nu = astro::trustRegion(x0,6);
        array<double,6> nu_arr{};
        for (int j=0;j<6;j++) nu_arr[j] = nu[j];
        nuOut.push_back(nu_arr);
    }

    // ------------------------------------------------------------
    // LOOP NODES 1..N-1
    // ------------------------------------------------------------
    for (int i=1;i<N;i++) {

        double t_i   = et0 + (i-1)*dt;
        double t_ip1 = et0 + (i)*dt;
        x = RK78Sc(6, x0, u, t_i, t_ip1, 1.0, Lsc,0, ap.gravOrd, dyn);
        if (flagProp == 1) {
            for (int j=0;j<6;j++)
                x0[j] = cons(x[j]) + DA(j+1);
        }
        else {
            // Use X0 columns directly (refinement propagation)
            for (int j=0;j<6;j++)
                x0[j] = X0.at(j,i) + DA(j+1);
        }
        // Save state & maps
        array<double,6> s{};
        for (int j=0;j<6;j++) s[j] = cons(x0[j]);
        state.push_back(s);

        auto STM = astro::stmDace(x,6,6);
        dynMaps.push_back(to_vec2m(STM,6,6));

        auto LLA = astro::llaMaps(x0, t_ip1, Lsc);
        auto LLAmap = astro::stmDace(LLA,6,3);

        vector<vector<double>> ll(2, vector<double>(3));
        for (int r=0;r<2;r++)
            for (int c=0;c<3;c++)
                ll[r][c] = LLAmap.at(r,c);
        llMaps.push_back(ll);

        llConst.push_back({ cons(LLA[0]), cons(LLA[1]) });

        auto nu = astro::trustRegion(x,6);
        array<double,6> nu_arr{};
        for (int j=0;j<6;j++) nu_arr[j] = nu[j];
        nuOut.push_back(nu_arr);
    }

    // ------------------------------------------------------------
    // WRITE JSON OUTPUT
    // ------------------------------------------------------------
    json jout;
    jout["state"]    = state;      // N × 6
    jout["dynMaps"]  = dynMaps;    // N × 6x6
    jout["llMaps"]   = llMaps;     // N × 2x3
    jout["llConst"]  = llConst;    // N × 2
    jout["nu"]       = nuOut;      // N × 6
    jout["timeMs"]   = now_ms() - t0;

    write_json("./data_sharing/geoOut.json", jout);

    clear_spice_kernels();
    return 0;

} catch (const std::exception& ex) {
    std::fprintf(stderr, "[findGeo] Error: %s\n", ex.what());
    return 1;
}