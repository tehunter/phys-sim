#ifndef PSYM_CORE_BOID_SIMPARAMETERS_H
#define PSYM_CORE_BOID_SIMPARAMETERS_H

namespace boid {

/** Simulation parameters for the boid simulation. */
struct SimParameters {

    // The distance from the wall at which point the wall starts to push back on the boid.
    static constexpr double WALL_MARGIN = 0.03;

    // The X position of the right wall (it's negative is the left wall).
    static constexpr double WALL_X = 0.8;
    // The Y position of the top wall (it's negative is the bottom wall).
    static constexpr double WALL_Y = 0.5;

    // The core grid parameters (how much the grid extends past the wall and how many cells in X/Y directions).
    static constexpr double GRID_PADDING = 0.20;
    static constexpr int NUM_GRID_CELLS = 50;

    // Computed constants for where the grid starts and ends.
    static constexpr double GRID_START_X = -WALL_X  - GRID_PADDING;
    static constexpr double GRID_END_X = WALL_X + GRID_PADDING;
    static constexpr double GRID_START_Y = -WALL_Y - GRID_PADDING;
    static constexpr double GRID_END_Y = WALL_Y + GRID_PADDING;

    // The per-cell widths and heights for the grid.
    static constexpr double GRID_CELL_WIDTH = (GRID_END_X - GRID_START_X) / NUM_GRID_CELLS;
    static constexpr double GRID_CELL_HEIGHT = (GRID_END_Y - GRID_START_Y) / NUM_GRID_CELLS;

    // Clicking within this distance of an existing goal removes it.
    static constexpr double GOAL_REMOVE_DIST = 0.03;

    /** The possible colors of boids. */
    enum Color { RED = 0, GREEN = 1, YELLOW = 2, BLUE = 3, LAST_COLOR_ELEMENT };

    /** The number of colors for boids. */
    static constexpr int NUM_COLORS = Color::LAST_COLOR_ELEMENT;

    /** When clicking, what should be spawned? **/
    enum ClickMode { CM_BOID1, CM_BOID5, CM_BOID20, CM_GOAL, CM_OBSTACLE, CM_PREDATOR };

    /** How should boids interact with the mouse? */
    enum MouseMode { MM_NONE, MM_AVOID, MM_FOLLOW };

    // The currently active color (see the enum).
    Color color;
    // The current click mode (see the enum).
    ClickMode click_mode;
    // The current mouse mode (see the enum).
    MouseMode mouse_mode;
    // If true, boids automatically avoid the cursor as if it were a predator.
    bool avoid_cursor;
    // How far does a boid look to find nearby boids to follow?
    double view_radius;
    // The base speed at which a boid moves.
    double speed;

    // If true, then boids will stay within an arena. Highly recommended to be enabled!
    bool walls_enabled;
    // How strongly walls push boids away from them.
    double wall_strength;
    // If true, then boids actively try to avoid being too close to other boids.
    bool separation_enabled;
    // The strength with which the boids react to the mouse.
    double mouse_strength;
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
    // If true, boids have an extrinistic desire to head towards goal locations.
    bool goals_enabled;
    // How strongly boids want to head towards a goal.
    double goal_strength;
    // If true, boids run away from predators with their dear lives.
    bool avoid_predators_enabled;
    // How strongly boids try to avoid predators.
    double avoid_predators_strength;
    // If true, boids avoid obstacles.
    bool avoid_obstacles_enabled;
    // How strongly boids avoid obstacles and predators.
    double avoid_obstacles_strength;

    // Predator view range (may be larger than a boid).
    double predator_view_radius;
    // Predator speed (may be faster than a boid).
    double predator_speed;
    // How strongly predators decide to chase nearby boids.
    double predator_chase_strength;

    SimParameters() {
        color = RED;
        click_mode = CM_BOID1;
        mouse_mode = MM_NONE;
        avoid_cursor = true;
        view_radius = 0.08;
        speed = 2.5e-4;
        walls_enabled = true;
        separation_enabled = true;
        follow_com_enabled = true;
        follow_vel_enabled = true;
        goals_enabled = true;
        avoid_obstacles_enabled = true;
        avoid_predators_enabled = true;

        wall_strength = 0.005;
        follow_com_strength = 0.005;
        follow_vel_strength = 0.005;
        separation_strength = 0.04;
        separation_threshold = 0.03;
        goal_strength = 0.01;
        avoid_obstacles_strength = 0.02;
        avoid_predators_strength = 0.05;

        predator_view_radius = 2 * view_radius;
        predator_speed = 1.6 * speed;
        predator_chase_strength = 0.05;
    }
};

}

#endif