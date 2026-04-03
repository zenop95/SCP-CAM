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

using namespace DACE;

int main() try {

    // --- DA init ---
    DA::init(2, 9);
    DA::setEps(1e-15);
    std::cerr << "[TEST] DA init OK" << std::endl;

    // --- SPICE ---
    furnsh_c("./data/kernel.txt");
    std::cerr << "[TEST] SPICE loaded OK" << std::endl;

    // --- AIDA flags: change these to isolate the problem ---
    int flags[3] = {2, 2, 2};  // flag1=atmo, flag2=srp, flag3=thirdbody
    // int flags[3] = {0, 0, 0};  // flag1=atmo, flag2=srp, flag3=thirdbody
    // Try combinations:
    //   {0,0,0} -> keplerian (should work)
    //   {1,0,0} -> atmosphere only
    //   {0,1,0} -> SRP only
    //   {0,0,1} -> third body only
    //   {1,1,1} -> full model

    const double Bfactor = 0.0044;
    const double SRPC    = 0.00262;
    const int    gravOrd = 5;
    const double Lsc     = 7000.0;  // km

    std::cerr << "[TEST] Building AIDAScaledDynamics with flags={"
              << flags[0] << "," << flags[1] << "," << flags[2]
              << "} gravOrd=" << gravOrd << std::endl;

    AIDAScaledDynamics<DA> dyn(
        "./data/gravmodels/egm2008",
        gravOrd, 
        flags,
        Bfactor,   
        SRPC       
    );
    std::cerr << "[TEST] AIDAScaledDynamics built OK" << std::endl;

    // --- Initial state: LEO orbit (scaled) ---
    // Position ~7000 km, velocity ~7.5 km/s, scaled by Lsc
    const double mu  = 398600.4418;
    const double Vsc = std::sqrt(mu / Lsc);
    AlgebraicVector<DA> x0(6);
    x0[0] = 1.0      + DA(1);   // x  [Lsc]
    x0[1] = 0.0      + DA(2);   // y
    x0[2] = 0.0      + DA(3);   // z
    x0[3] = 0.0      + DA(4);   // vx [Vsc]
    x0[4] = 1.0      + DA(5);   // vy  (circular orbit approx)
    x0[5] = 0.0      + DA(6);   // vz

    AlgebraicVector<DA> u(3, 0.0);  // zero control

    // --- Epoch (J2000 + 1 day in seconds) ---
    double et = 495943.194514429;

    // --- Single short propagation step ---
    double dt =  5828.51663768602*1;  // seconds (unscaled)
    double Tsc = Lsc / Vsc;
    double t0  = 0.0;
    double t1  = dt / Tsc;  // scaled time

    bool isNotKep = (gravOrd > 0) || (flags[0] > 0) || (flags[1] > 0) || (flags[2] > 0);
    std::cerr << "[TEST] isNotKep=" << isNotKep << std::endl;
    std::cerr << "[TEST] Starting propagation from t=" << t0 << " to t=" << t0 + dt << std::endl;

    AlgebraicVector<DA> x1 = RK78Sc(6, x0, u, et + t0, et + t1, 1.0, Lsc, false, isNotKep, dyn);

    std::cerr << "[TEST] Propagation OK" << std::endl;
    std::cerr << "[TEST] Final state (cons):" << std::endl;
    for (int i = 0; i < 6; ++i)
        std::cerr << "  x1[" << i << "] = " << cons(x1[i]) << std::endl;

    kclear_c();
    return 0;

} catch (const std::exception& ex) {
    std::cerr << "[TEST] Error: " << ex.what() << std::endl;
    return 1;
}