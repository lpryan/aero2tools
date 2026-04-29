core/
│
├── math/
│   ├── root_finding.cpp
│   ├── optimization.cpp
│
├── physics/
│   ├── isentropic.cpp
│   ├── oblique.cpp
│   ├── normal.cpp
│   ├── rayleigh.cpp
│   ├── nozzle.cpp
│
├── flow/
│   ├── State.hpp
│   ├── Flow.hpp
│   ├── Stage.hpp
│
├── stages/
│   ├── ObliqueStage.cpp
│   ├── ExpansionStage.cpp
│   ├── RayleighStage.cpp
│   ├── NozzleStage.cpp
│
├── solver/
│   ├── Solver.cpp
│   ├── Constraints.cpp
│
├── bindings/
│   ├── pybind11_interface.cpp

python/
│
├── aeropackage/
│   ├── __init__.py
│   ├── flow.py          # user-facing Flow wrapper
│   ├── stages.py        # enums / convenience
│   ├── state.py
│   ├── solver.py
│
├── examples/
│   ├── oblique.py
│   ├── nozzle.py
│   ├── rayleigh.py

tests/
├── homework/
│   ├── hw-1.py
│   ├── hw-2.py
│   ├── hw-3.py
│   ├── ...
│
├── exam/
│   ├── exam-1.py
│   ├── exam-2.py
│   ├── exam-3.py
│
├── feature/
│   ├── test_isentropic.cpp
│   ├── test_oblique.cpp
│   ├── test_flow_chain.cpp

.gitignore
makefile
pyproject.toml
README.md
setup.py