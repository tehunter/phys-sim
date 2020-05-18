#ifndef PSYM_CORE_BOID_SCENEOBJECTS_H
#define PSYM_CORE_BOID_SCENEOBJECTS_H

#include "SimParameters.h"

#include <Eigen/Core>

namespace boid {

/** A boid in the simulation. */
struct Boid {
    // The color of this boid.
    SimParameters::Color color;
    // The position of the boid.
    Eigen::Vector2d position;
    // The current velocity of the boid (this also provides it's orientation).
    Eigen::Vector2d velocity;
    // If true, this boid was marked for deletion.
    bool dead;

    inline Boid(SimParameters::Color color, Eigen::Vector2d pos, Eigen::Vector2d vel) : color(color), position(pos), velocity(vel), dead(false) {}
    inline Boid(SimParameters::Color color, Eigen::Vector2d pos) : Boid(color, pos, Eigen::Vector2d::Zero()) {}
};

struct Predator {
    static constexpr double BASE_SIZE_FACTOR = 1.6;

    // The position of this predator.
    Eigen::Vector2d position;
    // The velocity of this predator.
    Eigen::Vector2d velocity;

    // How many boids we've eaten.
    int eaten;

    inline Predator(Eigen::Vector2d pos, Eigen::Vector2d vel) : position(pos), velocity(vel), eaten(0) {}

    inline double size_factor() const {
        return BASE_SIZE_FACTOR + 0.05 * std::sqrt((double) eaten);
    }
};

/** A rectangular obstacle which has been placed in the simulation, and which boids will attempt to path around. */
struct Obstacle {
    // The position of the lower left corner of the obstacle.
    Eigen::Vector2d position;
    // The radius of the obstacle.
    double radius;

    inline Obstacle(Eigen::Vector2d pos, double radius) : position(pos), radius(radius) { }
};

}

#endif