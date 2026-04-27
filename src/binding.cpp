#include <pybind11/pybind11.h>

#include "add.hpp"

namespace py = pybind11;

PYBIND11_MODULE(aero2tools_core, m) {
    m.doc() = "C++ core for aero2tools";

    m.def("add", &cxx_add, "add two doubles");
}