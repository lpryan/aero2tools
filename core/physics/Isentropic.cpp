#include "State.hpp"
#include <cmath>

namespace Isentropic {

    double T0_T (State* state1) {

        double T0T = 1 + std::pow(state1->mach, 2) * CONFIG::GAMMA;
        
        
        // state1->mach






    }

}