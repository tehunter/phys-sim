#include "BoidHook.h"
#include <core/boid/BoidCore.h>

using namespace Eigen;

namespace boid {

BoidHook::BoidHook(BoidCore* core) : core_(core), PhysicsHook(core), params_(*core->getPointerToSimParameters()) {
}

BoidHook::~BoidHook() { }

void BoidHook::mouseClicked(double x, double y, int button) {
    message_mutex.lock();
    {
        MouseClick mc;
        mc.x = x;
        mc.y = y;
        mc.mode = params_.click_mode;
        mc.color = params_.color;
        mouseClicks_.push_back(mc);
    }
    message_mutex.unlock();
}

void BoidHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu) {
    if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Combo("Color", (int*)&params_.color, "Red\0Green\0Yellow\0Blue\0\0");
        ImGui::Combo("Click Mode", (int *)&params_.click_mode, "Spawn 1 Boid\0Spawn 5 Boids\0Spawn 20 Boids\0Spawner\0Goal\0\0");
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen)) {

    }
}

void BoidHook::tick() {
    enum ClickMode { CM_BOID1, CM_BOID5, CM_BOID20, CM_SPAWNER, CM_GOAL };
    message_mutex.lock();
    {
        while (!mouseClicks_.empty()) {
            MouseClick mc = mouseClicks_.front();
            mouseClicks_.pop_front();
            switch (mc.mode) {
                case SimParameters::ClickMode::CM_BOID1:
                    core_->spawn_boids(mc.x, mc.y, 1, mc.color);
                    break;
                case SimParameters::ClickMode::CM_BOID5:
                    core_->spawn_boids(mc.x, mc.y, 5, mc.color);
                    break;
                case SimParameters::ClickMode::CM_BOID20:
                    core_->spawn_boids(mc.x, mc.y, 20, mc.color);
                    break;
                case SimParameters::ClickMode::CM_SPAWNER:
                    core_->create_spawner(mc.x, mc.y, mc.color);
                    break;
                case SimParameters::ClickMode::CM_GOAL:
                    core_->create_goal(mc.x, mc.y, mc.color);
                    break;
            }
        }
    }
    message_mutex.unlock();
}

}