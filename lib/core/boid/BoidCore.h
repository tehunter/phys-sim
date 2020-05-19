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
#include <random>

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

    // Spawn a predator at the given location.
    void spawn_predator(double x, double y);

    // Create a boid goal at the given location with the given color.
    void create_goal(double x, double y, SimParameters::Color color);

    // Create an obstacle with the given radius.
    void create_obstacle(double x, double y);

    // Update the current mouse position (for use in follow/avoid behaviors).
    void update_mouse_position(double x, double y);

    // Obtain a random initial velocity for a newly spawned boid.
    Eigen::Vector2d random_initial_velocity();

    // Obtain the cell corresponding to the given map location.
    Eigen::Vector2i grid_cell(Eigen::Vector2d pos) const;

    // Find all of the neighbors to the boid with the given ID, placing them in the output vector.
    int fill_neighbors(Eigen::Vector2d position, double radius, int ignore, std::vector<int>& output) const;

    // Obtain the goal location for the given boid color, if present (the first boolean determines if there is a goal to
    // follow).
    std::pair<bool, Eigen::Vector2d> goal(SimParameters::Color color) const;

    // Obtains the velocity towards the center of mass of the given neighbors.
    Eigen::Vector2d velocity_com(const Boid& us, std::vector<int>& neighbors) const;
    // Obtains the velocity which matches the velocity of the given neighbors.
    Eigen::Vector2d velocity_vel(const Boid& us, std::vector<int>& neighbors) const;
    // Obtains the velocity which avoids nearby neighbors which are too close.
    Eigen::Vector2d velocity_avoid(const Boid& us, std::vector<int>& neighbors) const;
    // Obtains the velocity for avoiding running into walls (and off screen).
    Eigen::Vector2d velocity_wall(Eigen::Vector2d position) const;
    // Obtains the velocity for heading towards a goal.
    Eigen::Vector2d velocity_goal(const Boid& us) const;
    // Obtains the velocity for avoiding nearby predators.
    Eigen::Vector2d velocity_avoid_preds(const Boid& us) const;
    // Obtains the velocity for avoiding obstacles.
    Eigen::Vector2d velocity_obstacle(Eigen::Vector2d position) const;

    // Internal method which returns a mesh representing the current state.
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> getCurrentMesh() const override;

    // Utility methods used by the visualizer to obtain the simulation parameters.
    inline std::shared_ptr<SimParameters> getPointerToSimParameters() { return params_; }
private:
    // The list of all boids.
    std::vector<Boid> boids;

    // The list of all predators.
    std::vector<Predator> predators;

    // A basic spatial partitioning data strcture which is a 2D grid.
    std::vector<std::vector<std::vector<int>>> grid;

    // The current location of the mouse. Used for follow/avoid mechanisms.
    Eigen::Vector2d mouse_location;

    // The goal for each boid; the boolean determines if a goal has been set, and the position if the goals location.
    std::pair<bool, Eigen::Vector2d> goals[SimParameters::NUM_COLORS];

    // The list of obstacles which boids must attempt to path around.
    std::vector<Obstacle> obstacles;

    // Random number generator.
    std::mt19937 rng;

    // The parameters which control the simulation.
    std::shared_ptr<SimParameters> params_;
};

}

#endif