#include <nlohmann/json.hpp>
#include <vector>
#include <array>
#include <fstream>
#include <stdexcept>
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"
#include "dynorb/AIDA.h"
#include "dynorb/AIDAwrappers.h"
#include "prop_utils.h"

using json = nlohmann::json;
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
    // ------------------ Timing ------------------
    long t0 = now_ms();
    long time_acc_ms = 0;

    // ------------------ Read JSON input ------------------
    auto jin = read_json("./data_sharing/initial_state.json");

    const int N = jin.at("N").get<int>();
    const std::vector<double> t = jin.at("t").get<std::vector<double>>();
    const double et   = jin.at("et").get<double>();
    const double aMax = jin.at("uMax").get<double>();
    const double Lsc  = jin.at("scaling").get<double>();
    const int order   = jin.at("order").get<int>();
    const int flagRtn = jin.at("flagRtn").get<int>();
    const int flagProp= jin.at("flagProp").get<int>();

    // Control and trajectory grids
    std::vector<std::vector<double>> xdum;          // N x 6 (only used if flagProp==2)
    const std::vector<std::vector<double>> udum = jin.at("u").get<std::vector<std::vector<double>>>();

    // ------------------ DA / AIDA / SPICE init ------------------
    DA::init(order, 9);
    DA::setEps(1e-15);

    const auto aida = read_aida_params_json("./data_sharing/AIDA_init.json");
    load_spice_kernel("./data/kernel.txt");
    auto coeff = compute_aida_coeff(aida);
    int flags[3] = { aida.flag1, aida.flag2, aida.flag3 };
    AIDAScaledDynamics<DA> dyn(
        "./data/gravmodels/egm2008",
        aida.gravOrd,
        flags,
        coeff.Bfactor,
        coeff.SRPC
    );
    bool isNotKep = (aida.gravOrd > 1) || (aida.flag1 > 0) || (aida.flag2 > 0) || (aida.flag3 > 0);
    // ------------------ State / control symbols ------------------
    AlgebraicVector<DA> x0(6), x(6), u(3), lla(3);

    // Initialize x0 for node 0
    if (flagProp == 1) {
        // Read x0 from JSON and make it DA
        const auto x0_dum = jin.at("x0").get<std::vector<double>>();
        for (int j = 0; j < 6; ++j) {
            x0[j] = x0_dum[j] + DA(1)*0;
        }
        if (t[0] != 0){
            x0 = RK78Sc(6, x0, u, et, et + t[0], 1.0, Lsc, 0, isNotKep, dyn);
        }
    } else {
        // Read full trajectory and take node 0 as initial DA
        xdum = jin.at("oldTraj").get<std::vector<std::vector<double>>>();
        for (int j = 0; j < 6; ++j) {
            x0[j] = xdum[0][j] + DA(1)*0;
        }
    }
    for (int j = 0; j < 6; ++j) {
        x0[j] = cons(x0[j]) + DA(j+1);
    }
    // ------------------ Output containers ------------------
    std::vector<std::array<double,6>> state;                     // N × (6)
    std::vector<std::vector<std::vector<double>>> dynMaps;       // N × (6×9)
    std::vector<std::array<double,3>> latConst;                  // N × (3)
    std::vector<std::vector<std::vector<double>>> latMaps;       // N × (3×9)
    std::vector<std::array<double,9>> trusts;                    // N × (9) only if order>1

    // ------------------ Node 0 outputs ------------------
    {
        // lla maps and STMs at node 0
        lla          = astro::llaMaps(x0, et, Lsc);
        auto STM     = astro::stmDace(x0, 9, 6);   // returns 6x9
        auto llaSTM  = astro::stmDace(lla, 9, 3);  // returns 3x9

        // Constant state at node 0
        std::array<double,6> s{};
        for (int j = 0; j < 6; ++j) s[j] = cons(x0[j]);
        state.push_back(s);

        // Store matrices with correct shapes
        dynMaps.push_back(to_vec2m(STM,    6, 9));
        latConst.push_back({ cons(lla[0]), cons(lla[1]), cons(lla[2]) });
        latMaps.push_back(to_vec2m(llaSTM, 3, 9));

        if (order > 1) {
            auto nu = astro::trustRegion(x0, 9);
            std::array<double,9> a{};
            for (int k = 0; k < 9; ++k) a[k] = nu[k];
            trusts.push_back(a);
        }
    }

    long t1 = now_ms();
    time_acc_ms += now_ms() - t1;
    // ------------------ Nodes 1..N-1 ------------------
    for (int i = 1; i < N; ++i) {
        // Propagate from node i-1 to i using previous x0 and control u
        // Control for node 0
        for (int j = 0; j < 3; ++j) {
            u[j] = udum[i-1][j] + DA(j+7);
        }

        x = RK78Sc(6, x0, u, et + t[i-1], et + t[i], aMax, Lsc, flagRtn, isNotKep, dyn);
        // Outputs for node i
        auto STM    = astro::stmDace(x,   9, 6);   // 6x9
        auto llaSTM = astro::stmDace(lla, 9, 3);   // 3x9
        
        if (flagProp == 1) {
            // Update x0 from propagated x (freeze then add DA variables)
            for (int j = 0; j < 6; ++j) x0[j] = cons(x[j]) + DA(j+1);
        } else {
            // Use provided trajectory for node i, plus DA variables
            for (int j = 0; j < 6; ++j) x0[j] = xdum[i][j] + DA(j+1);
        }
        lla         = astro::llaMaps(x0, et + t[i], Lsc);

        std::array<double,6> s{};
        for (int j = 0; j < 6; ++j) s[j] = cons(x[j]);
        state.push_back(s);

        dynMaps.push_back(to_vec2m(STM,    6, 9));
        latConst.push_back({ cons(lla[0]), cons(lla[1]), cons(lla[2]) });
        latMaps.push_back(to_vec2m(llaSTM, 3, 9));

        if (order > 1) {
            auto nu = astro::trustRegion(x, 9);
            std::array<double,9> a{};
            for (int k = 0; k < 9; ++k) a[k] = nu[k];
            trusts.push_back(a);
        }
    }

    // ------------------ Write output JSON ------------------
    time_acc_ms += now_ms() - t0;

    json jout;
    jout["state"]        = state;       // 6 per node
    jout["dynMaps"]      = dynMaps;     // each 6x9
    jout["latlonConst"]  = latConst;    // each 3
    jout["latlonMaps"]   = latMaps;     // each 3x9
    if (order > 1) jout["trustRegion"] = trusts;
    jout["timeMs"]       = time_acc_ms;

    write_json("./data_sharing/out_prop.json", jout);
    clear_spice_kernels();
    return 0;

} catch (const std::exception& ex) {
    std::fprintf(stderr, "[propAida] Error: %s\n", ex.what());
    return 1;
}