#ifndef PSYM_CORE_BOID_BOIDCORE_H
#define PSYM_CORE_BOID_BOIDCORE_H

#include "../PhysicsCore.h"
#include "SimParameters.h"

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <memory>
#include <stdint.h>

namespace boid {

class BoidCore : public PhysicsCore {
public:
    BoidCore();
    ~BoidCore();

    // Initialize an empty simulation.
    virtual void initSimulation() override;

    // Simulate the simulation by a single discrete step.
    virtual bool simulateOneStep() override;


    // Spawn the given number of boids of the given color at the given location.
    void spawn_boids(double x, double y, int count, SimParameters::Color color);

    // Create a boid spawner at the given location with the given color.
    void create_spawner(double x, double y, SimParameters::Color color);

    // Create a boid goal at the given location with the given color.
    void create_goal(double x, double y, SimParameters::Color color);


    // Internal method which returns a mesh representing the current state.
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> getCurrentMesh() const override;

    // Utility methods used by the visualizer to obtain the simulation parameters.
    inline std::shared_ptr<SimParameters> getPointerToSimParameters() { return params_; }
private:
    std::shared_ptr<SimParameters> params_;
};

}

#endif