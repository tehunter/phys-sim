#include "../IglVisualizer.h"
#include <thread>
#include <memory>

static PhysicsHook* hook = nullptr;

namespace boid {

class BoidCore;
class BoidHook;

class BoidVisualizer : public IglVisualizer {
    std::unique_ptr<BoidCore> core_;
    std::unique_ptr<BoidHook> hook_;
public:
    BoidVisualizer();
    ~BoidVisualizer();
};

}