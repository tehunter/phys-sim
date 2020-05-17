#ifndef PSYM_CORE_BOID_BOIDCORE_H
#define PSYM_CORE_BOID_BOIDCORE_H

#include "../PhysicsCore.h"
#include "SimParameters.h"
#include "SceneObjects.h"

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <memory>
#include <vector>
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
    // Per-color lists of boids.
    std::vector<Boid> boids[SimParameters::NUM_COLORS];

    // The goal for each boid; the boolean determines if a goal has been set, and the position if the goals location.
    std::pair<bool, Eigen::Vector2d> goals[SimParameters::NUM_COLORS];

    // The list of obstacles which boids must attempt to path around.
    std::vector<Obstacle> obstacles;

    std::shared_ptr<SimParameters> params_;
};

}

#endif