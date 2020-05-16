#include "BoidHook.h"
#include "BoidVisualizer.h"
#include <core/boid/BoidCore.h>

namespace boid {

BoidVisualizer::BoidVisualizer() {
    core_.reset(new BoidCore);
    hook_.reset(new BoidHook(core_.get()));
    hook_->reset();
    init(core_.get(), hook_.get());
}

BoidVisualizer::~BoidVisualizer() { }

}