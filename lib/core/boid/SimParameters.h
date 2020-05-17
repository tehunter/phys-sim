#ifndef PSYM_CORE_BOID_SIMPARAMETERS_H
#define PSYM_CORE_BOID_SIMPARAMETERS_H

namespace boid {

/** Simulation parameters for the boid simulation. */
struct SimParameters {

    /** The possible colors of boids. */
    enum Color { RED = 0, GREEN = 1, YELLOW = 2, BLUE = 3, LAST_COLOR_ELEMENT };

    /** The number of colors for boids. */
    static constexpr int NUM_COLORS = Color::LAST_COLOR_ELEMENT;

    /** When clicking, what should be spawned? **/
    enum ClickMode { CM_BOID1, CM_BOID5, CM_BOID20, CM_SPAWNER, CM_GOAL };

    // The currently active color (see the enum).
    Color color;
    // The current click mode (see the enum).
    ClickMode click_mode;
    // If true, boids automatically avoid the cursor as if it were a predator.
    bool avoid_cursor;
    // The height and width of the arena the boids live in; they will actively try to avoid walls.
    double arena_width, arena_height;
    // How far does a boid look to find nearby boids to follow?
    double view_radius;
    // The base speed at which a boid moves.
    double speed;

    // If true, then boids will stay within an arena. Highly recommended to be enabled!
    bool arena_enabled;
    // If true, then boids actively try to avoid being too close to other boids.
    bool separation_enabled;
    // How strongly to weight the separation force.
    double separation_strength;
    // At what distance does separation start to kick in? Boids will tend to space themself out by this threshold.
    double separation_threshold;
    // If true, boids tend to move towards the center of mass of surrounding boids.
    bool follow_com_enabled;
    // How strongly to weight the following center of mass force.
    double follow_com_strength;
    // If true, boids tend to move in the same direction as other boids are currently moving.
    bool follow_vel_enabled;
    // How strongly to weight the following velocity force.
    double follow_vel_strength;

    // If true, boids have inertia in their movement which slows down turns and reactions.
    bool inertia_enabled;
    // How strongly the boid wants to continue going in it's current direction.
    double inertia_strength;

    SimParameters() {
        color = RED;
        click_mode = CM_BOID1;
        avoid_cursor = true;
        arena_width = arena_height = 1.0;
        view_radius = 0.15;
        speed = 4e-5;
        arena_enabled = true;
        separation_enabled = true;
        follow_com_enabled = true;
        follow_vel_enabled = true;
        inertia_enabled = true;

        follow_com_strength = 1e-2;
        follow_vel_strength = 1e-1;
        inertia_strength = 1.0;
        separation_strength = 1.0;
        separation_threshold = 0.08;
    }
};

}

#endif