#ifndef STATE_H
#define STATE_H

#include "../physics/Config.hpp"

struct State {

    double mach = 1.0;
    
    double T = 0.0;
    double T0 = 0.0;
    double P = 0.0;
    double P0 = 0.0;

    State(double M) : mach(M) {}

    State(double M, double T, double P) : mach(M), T(T), P(P) {}

};

#endif//STATE_H