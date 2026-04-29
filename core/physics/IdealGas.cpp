#include "../physics/IdealGas.hpp"
#include "../physics/Config.hpp"

namespace IdealGas {

    double pres (double temp, double dens) {
        return dens * CONFIG::R * temp;
    }

    double temp (double pres, double dens) {
        
    }



    double dens (double temp, double pres);

}