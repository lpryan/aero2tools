#ifndef CONFIG_H
#define CONFIG_H

namespace CONFIG {

    double GAMMA = 1.4;
    double R = 287; // J/(kg*K)

    double CP() { return GAMMA * R / (GAMMA - 1); }

    double CV() { return R / (GAMMA - 1); }

};

#endif//CONFIG_H