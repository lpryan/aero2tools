CXX = g++
CXXFLAGS = -O3 -Wall -shared -std=c++17 -fPIC -Iinclude

PYTHON_INCLUDE := $(shell python3 -m pybind11 --includes)
PYTHON_EXT_SUFFIX := $(shell python3-config --extension-suffix)

SRC = src
PKG = aero2tools/cpp

TARGET = $(PKG)/aero2tools_core$(PYTHON_EXT_SUFFIX)

MATH_LIB = $(PKG)/libcxx_math.so
OPTIMIZE_LIB = $(PKG)/libcxx_optimize.so

all: $(TARGET)

$(OPTIMIZE_LIB): $(SRC)/optimize/optimize.cpp
	mkdir -p $(PKG)
	$(CXX) $(CXXFLAGS) -shared $^ -o $@

$(MATH_LIB): $(SRC)/add.cpp
	mkdir -p $(PKG)
	$(CXX) $(CXXFLAGS) -shared $^ -o $@


$(TARGET): $(SRC)/binding.cpp $(OPTIMIZE_LIB) $(MATH_LIB)
	mkdir -p $(PKG)
	$(CXX) $(CXXFLAGS) $(SRC)/binding.cpp -o $@ \
	    $(PYTHON_INCLUDE) -L$(PKG) -lcxx_math -lcxx_optimize \
	    -Wl,-rpath,'$$ORIGIN'

clean:
	rm -f $(PKG)/*.so