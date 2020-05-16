#include "boid.h"
#include <core/boid/BoidCore.h>

#if PSIM_ENABLE_VISUALIZER
#include <vis/boid/BoidVisualizer.h>
#endif

namespace boid {

void define_module(py::module& m) {
#if PSIM_ENABLE_VISUALIZER
	py::class_<BoidVisualizer, IglVisualizer>(m, "BoidVisualizer", py::module_local())
		.def(py::init<>());
#endif
}

};
