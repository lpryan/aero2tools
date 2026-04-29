#ifndef IDEALGAS_H
#define IDEALGAS_H

namespace IdealGas {

    double pres (double temp, double dens);

    double temp (double pres, double dens);

    double dens (double temp, double pres);

};

#endif//IDEALGAS_H