#include "prop_utils.h"
#include <nlohmann/json.hpp>
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"

using namespace props::io;
using namespace DACE;
using json = nlohmann::json;

int main() try {

    auto in = read_input_tca_json("./data_sharing/initial_state.json");
    long t0 = now_ms();

    DA::init(1,7);
    DA::setEps(1e-15);

    AlgebraicVector<DA> xp0(6), xs0(6), xpf(6), xsf(6), xrel(6);
    DA dt, tca;

    std::vector<double> tca_out;
    tca_out.reserve(in.M);

    for (int j=0; j<in.M; j++) {

        dt = 0.0 + DA(7);

        for (int i=0;i<6;i++) {
            xp0[i] = in.xdump.at(i,j);
            xs0[i] = in.xdums.at(i,j) + DA(i+1);
        }

        xpf = KeplerProp(xp0,dt,1.0);
        xsf = KeplerProp(xs0,dt,1.0);

        xrel = xpf - xsf;

        tca = astro::findTCA(xrel,7);
        tca_out.push_back( cons(tca) );
    }

    json jout;
    jout["tca"]    = tca_out;
    jout["timeMs"] = now_ms() - t0;

    write_json("./data_sharing/tcaOut.json", jout);
    return 0;

} catch(const std::exception& ex) {
    std::fprintf(stderr,"[TCA] Error: %s\n", ex.what());
    return 1;
}