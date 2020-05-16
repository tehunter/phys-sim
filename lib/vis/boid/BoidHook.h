#ifndef PSYM_VIS_BOID_BOIDHOOK_H
#define PSYM_VIS_BOID_BOIDHOOK_H

#include "../PhysicsHook.h"
#include <core/boid/SimParameters.h>
#include <deque>

namespace boid {

class BoidCore;

struct MouseClick {
    double x;
    double y;
    SimParameters::ClickMode mode;
    SimParameters::Color color;
};

class BoidHook : public PhysicsHook {
public:
    BoidHook(BoidCore* core);
    ~BoidHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu) override;

    virtual void mouseClicked(double x, double y, int button) override;

    virtual void tick();

private:
    BoidCore* core_;
    SimParameters& params_;
    double time_;

    std::mutex message_mutex;
    std::deque<MouseClick> mouseClicks_;
};

}

#endif
