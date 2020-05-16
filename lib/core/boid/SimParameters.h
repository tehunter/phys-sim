#ifndef PSYM_CORE_BOID_SIMPARAMETERS_H
#define PSYM_CORE_BOID_SIMPARAMETERS_H

namespace boid {

/** Simulation parameters for the boid simulation. */
struct SimParameters {

    /** The possible colors of boids. */
    enum Color { RED = 0, GREEN = 1, YELLOW = 2, BLUE = 3, NUM_COLORS };

    /** When clicking, what should be spawned? **/
    enum ClickMode { CM_BOID1, CM_BOID5, CM_BOID20, CM_SPAWNER, CM_GOAL };

    // The currently active color (see the enum).
    Color color;
    // The current click mode (see the enum).
    ClickMode click_mode;
    // If true, boids automatically avoid the cursor as if it were a predator.
    bool avoid_cursor;

    SimParameters() {
        color = RED;
        click_mode = CM_BOID1;
        avoid_cursor = true;
    }
};

}

#endif