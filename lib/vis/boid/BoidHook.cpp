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

void BoidHook::mouseMoved(double x, double y) {
    message_mutex.lock();
    {
        MouseMove mm;
        mm.x = x;
        mm.y = y;
        mouse_moves.push_back(mm);
    }
    message_mutex.unlock();
}

void BoidHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu) {
    if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Combo("Color", (int*)&params_.color, "Red\0Green\0Yellow\0Blue\0\0");
        ImGui::Combo("Click Mode", (int *)&params_.click_mode, "Spawn 1 Boid\0Spawn 5 Boids\0Spawn 20 Boids\0Goal\0Obstacle\0Predator\0\0");
        ImGui::Combo("Mouse Mode", (int*) &params_.mouse_mode, "None\0Avoid\0Follow\0\0");
        ImGui::InputDouble("Obstacle Radius", &params_.obstacle_radius);
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputDouble("Speed", &params_.speed);
        ImGui::InputDouble("View Radius", &params_.view_radius);
    }
    if (ImGui::CollapsingHeader("Rules", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Walls Enabled", &params_.walls_enabled);
        ImGui::InputDouble("Wall Strength", &params_.wall_strength);
        ImGui::Checkbox("Separation Enabled", &params_.separation_enabled);
        ImGui::InputDouble("Separation Strength", &params_.separation_strength);
        ImGui::InputDouble("Separation Distance", &params_.separation_threshold);
        ImGui::Checkbox("COM Enabled", &params_.follow_com_enabled);
        ImGui::InputDouble("COM Strength", &params_.follow_com_strength);
        ImGui::Checkbox("VEL Enabled", &params_.follow_vel_enabled);
        ImGui::InputDouble("VEL Strength", &params_.follow_vel_strength);
        ImGui::Checkbox("Goals Enabled", &params_.goals_enabled);
        ImGui::InputDouble("Goal Strength", &params_.goal_strength);
        ImGui::Checkbox("Avoid Pred. Enabled", &params_.avoid_predators_enabled);
        ImGui::InputDouble("Avoid Pred. Strength", &params_.avoid_predators_strength);
        ImGui::Checkbox("Avoid Obs. Enabled", &params_.avoid_obstacles_enabled);
        ImGui::InputDouble("Avoid Obs. Strength", &params_.avoid_obstacles_strength);
    }
    if (ImGui::CollapsingHeader("Predators", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputDouble("Predator View Range", &params_.predator_view_radius);
        ImGui::InputDouble("Predator Speed", &params_.predator_speed);
        ImGui::InputDouble("Chase Strength", &params_.predator_chase_strength);
    }
}

void BoidHook::tick() {
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
                case SimParameters::ClickMode::CM_GOAL:
                    core_->create_goal(mc.x, mc.y, mc.color);
                    break;
                case SimParameters::ClickMode::CM_OBSTACLE:
                    core_->create_obstacle(mc.x, mc.y);
                    break;
                case SimParameters::ClickMode::CM_PREDATOR:
                    core_->spawn_predator(mc.x, mc.y);
                    break;
            }
        }

        while (!mouse_moves.empty()) {
            MouseMove mm = mouse_moves.front();
            mouse_moves.pop_front();
            core_->update_mouse_position(mm.x, mm.y);
        }
    }
    message_mutex.unlock();
}

}
