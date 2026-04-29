#ifndef SYSTEMVARIABLE_H
#define SYSTEMVARIABLE_H

#include <unordered_map>
#include <string>
#include <stdexcept>
#include <vector>

#include "Variable.hpp"

class SystemVariables {
private:
    std::unordered_map<std::string, Variable> vars;

public:

    void add(const std::string& name, double value, VarType type) {
        vars[name] = Variable(name, value, type);
    }

    double& operator[](const std::string& name) {
        if (!vars.count(name)) {
            throw std::runtime_error("Variable not found: " + name);
        }
        return vars[name].value;
    }

    const double& at(const std::string& name) const {
        return vars.at(name).value;
    }

    bool exists(const std::string& name) const {
        return vars.count(name);
    }

    std::vector<std::string> getActiveUnknowns() const {
        std::vector<std::string> keys;
        for (const auto& [k, v] : vars) {
            if (v.type == VarType::STATE && v.active) {
                keys.push_back(k);
            }
        }
        return keys;
    }
};

#endif // SYSTEMVARIABLE_H