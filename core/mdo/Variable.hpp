#ifndef VARIABLE_H
#define VARIABLE_H

#include <string>
#include <unordered_map>

enum class VarType {
    DESIGN,     // controlled by optimizer
    STATE,      // solved internally
    TARGET,     // constants
};

struct Variable {
    std::string name;
    double value = 0.0;
    VarType type = VarType::STATE;
    bool active = true;

    Variable() = default;

    Variable(std::string n, double v, VarType t) : name(std::move(n)), value(v), type(t) {}
};

#endif // VARIABLE_H