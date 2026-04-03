#include <dace/dace.h>
#include <nlohmann/json.hpp>
#include <vector>
#include <array>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include "prop_utils.h"

#include "astro/AstroRoutines.h"

using json = nlohmann::json;
using namespace DACE;
using namespace std;
using namespace std::chrono;
using namespace props::io;

int main() try {

    long t0 = now_ms();

    // --------- Read input JSON ---------
    auto jin = read_json("./data_sharing/ipcIn.json");
    int M    = jin.at("M").get<int>();
    int type = jin.at("type").get<int>();

    // Accept scalar or array for weights/R
    auto weights = read_vec(jin, "weights");
    auto R       = read_vec(jin, "R");

    // Accept [x,y,z] (M=1) or [[x,y,z], ...] (M>=1)
    auto r_in    = read_r_list(jin, "r", 1, M);

    // Accept single 3x3 OR list of 3x3
    auto P_in    = read_mat3_list(jin, "P", M);
    auto e2b_in  = read_mat3_list(jin, "e2b", M);

    // Basic size checks (after normalization these must match)
    if ((int)weights.size() != M) throw std::runtime_error("weights size != M");
    if ((int)R.size()       != M) throw std::runtime_error("R size != M");
    if ((int)r_in.size()    != M) throw std::runtime_error("r size != M");
    if ((int)P_in.size()    != M) throw std::runtime_error("P size != M");
    if ((int)e2b_in.size()  != M) throw std::runtime_error("e2b size != M");

    // --------- DA init ---------
    DA::init(2, 3*M);
    DA::setEps(1e-15);

    // Containers
    AlgebraicVector<DA> r(3), rB(3), rB2(2), poc(M);
    AlgebraicMatrix<double> cov(3,3), e2b(3,3), covB(3,3), covB2(2,2);

    // Compute poc[i], accumulate total
    DA noCollisions = 1.0;

    int w = 0;  // variable index mapping (3 per object)
    for (int i=0; i<M; ++i) {

        // r (3) with DA variables
        for (int j=0; j<3; ++j) r[j] = r_in[i][j] + DA(w + (j+1));
        w += 3;

        // cov 3x3
        for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k)
                cov.at(j,k) = P_in[i][j][k];

        // e2b 3x3
        for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k)
                e2b.at(j,k) = e2b_in[i][j][k];

        // rotate to B-plane
        covB = e2b * cov * e2b.transpose();
        rB   = e2b * r;

        // 2D reduction (0 and 2 components)
        rB2[0] = rB[0]; rB2[1] = rB[2];
        covB2.at(0,0)=covB.at(0,0); covB2.at(0,1)=covB.at(0,2);
        covB2.at(1,0)=covB.at(2,0); covB2.at(1,1)=covB.at(2,2);

        // per-object PoC
        if      (type==0) poc[i] = astro::ConstPoC (rB2, covB2, R[i]) * weights[i];
        else if (type==1) poc[i] = astro::MaxPoC   (rB2, covB2, R[i]) * weights[i];
        else if (type==2) poc[i] = astro::ChanPoC  (rB2, covB2, R[i], 2) * weights[i];
        else if (type==3) poc[i] = astro::AlfanoPoC(rB2, covB2, R[i]) * weights[i];
        else if (type==4) poc[i] = rB2[0]*rB2[0] + rB2[1]*rB2[1];  // miss distance
        else throw std::runtime_error("PoC type must be in [0..4]");

        noCollisions = noCollisions * (1 - poc[i]);
    }

    DA poc_tot = 1.0 - noCollisions;

    // Trust region on total PoC
    AlgebraicVector<double> nu = astro::trustRegionScalar(poc_tot, 3*M);

    // Build pcMaps (M x 3) — derivative mapping
    std::vector<std::array<double,3>> pcMaps(M);
    int varIndex = 0;
    for (int i=0; i<M; ++i) {
        for (int j=0; j<3; ++j) {
            double d;
            if (type == 4)
                d = DACE::cons( poc[i].deriv(varIndex + 1) );
            else
                d = DACE::cons( poc_tot.deriv(varIndex + 1) );
            pcMaps[i][j] = d;
            varIndex++;
        }
    }

    // Flatten trustRegion (3*M) to match MATLAB reshape(...,3,M)
    std::vector<double> trust_flat(3*M);
    for (int k=0; k<3*M; ++k) trust_flat[k] = nu[k];

    // Write JSON
    json jout;
    jout["pcMaps"]      = pcMaps;        // Mx3 (d/d scaled coords)
    jout["PoC"]         = DACE::cons(poc_tot);
    jout["trustRegion"] = trust_flat;    // length 3*M
    jout["timeMs"]      = now_ms() - t0;

    write_json("./data_sharing/pocOut.json", jout);
    return 0;

} catch (const std::exception& ex) {
    std::fprintf(stderr, "[pocMaps] Error: %s\n", ex.what());
    return 1;
}