#ifndef PSIM_LIB_PYTHON_BOID_H
#define PSIM_LIB_PYTHON_BOID_H

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

namespace boid {
void define_module(py::module& m);
};

#endif
