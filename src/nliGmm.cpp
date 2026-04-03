#include "prop_utils.h"
#include <nlohmann/json.hpp>
#include "astro/AstroLibrary.h"
#include "astro/AstroRoutines.h"

using namespace props::io;
using namespace DACE;
using json = nlohmann::json;

int main() try {

    auto in = read_input_trustregion_json("./data_sharing/nli.json");

    DA::init(2,6);
    DA::setEps(1e-15);

    AlgebraicVector<DA> x0(6), x(6);
    for (int j=0;j<6;j++)
        x0[j] = in.x0[j] + DA(j+1);

    x = KeplerProp(x0, in.dt, 1.0);

    auto nu = astro::trustRegion(x,6);

    json jout;
    jout["trustRegion"] =
       {nu[0], nu[1], nu[2], nu[3], nu[4], nu[5]};
    write_json("./data_sharing/trustRegion.json", jout);
    return 0;

} catch(const std::exception& ex) {
    std::fprintf(stderr,"[TrustRegion] Error: %s\n", ex.what());
    return 1;
}