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

    inline Boid(SimParameters::Color color, Eigen::Vector2d pos, Eigen::Vector2d vel) : color(color), position(pos), velocity(vel) {}
    inline Boid(SimParameters::Color color, Eigen::Vector2d pos) : Boid(color, pos, Eigen::Vector2d::Zero()) {}
};

/** A rectangular obstacle which has been placed in the simulation, and which boids will attempt to path around. */
struct Obstacle {
    // The position of the lower left corner of the obstacle.
    Eigen::Vector2d position;
    // The size (width, height) of the obstacle.
    Eigen::Vector2d size;

    inline Obstacle(Eigen::Vector2d pos, Eigen::Vector2d size) : position(pos), size(size) { }
};

/** A spawner which endlessly creates boids at a given rate. */
struct Spawner {
    // The color of boid that this spawner creates.
    SimParameters::Color color;
    // The location of this spawner.
    Eigen::Vector2d position;
    // How frequently (i.e., how many boids/second) this spawner creates boids.
    double rate;

    inline Spawner(SimParameters::Color color, Eigen::Vector2d pos, double rate) : color(color), position(pos), rate(rate) { }
};

}

#endif