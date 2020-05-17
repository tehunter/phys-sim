#include "BoidCore.h"
#include "SimParameters.h"

#include <iostream>

#include <Eigen/Geometry>

using namespace Eigen;

namespace boid {

// Trivial utility method which promotes a 2d vector to a 3d vector with zero-ish z-component.
Vector3d z0(const Vector2d input) {
    return Vector3d(input(0), input(1), 0);
}

BoidCore::BoidCore() {
    params_ = std::make_shared<SimParameters>();
}

BoidCore::~BoidCore() { }

void BoidCore::initSimulation() {
    for (int color = 0; color < SimParameters::NUM_COLORS; color++) {
        this->boids[color].clear();
        this->goals[color] = std::make_pair(false, Vector2d::Zero());
        this->obstacles.clear();
    }

    // Place down some initial boids.
    this->boids[0].emplace_back(SimParameters::RED, Vector2d(0.0, 0.0), Vector2d(1.0, -1.0));
    this->boids[0].emplace_back(SimParameters::RED, Vector2d(0.05, 0.0), Vector2d(1.0, -1.0));
    this->boids[0].emplace_back(SimParameters::RED, Vector2d(-0.05, 0.0), Vector2d(1.0, -1.0));
    this->boids[0].emplace_back(SimParameters::RED, Vector2d(0.0, 0.05), Vector2d(1.0, -1.0));
    this->boids[0].emplace_back(SimParameters::RED, Vector2d(0.0, -0.05), Vector2d(1.0, -1.0));
    this->boids[1].emplace_back(SimParameters::GREEN, Vector2d(0.5, 0.0), Vector2d(0.0, 1.0));
    this->boids[2].emplace_back(SimParameters::GREEN, Vector2d(0.5, 0.5), Vector2d(-1.0, 1.0));
}

bool BoidCore::simulateOneStep() {
    // TODO: Boid states are updated live (so one boid may see an updated velocity, while another may not);
    // this may result in some wierd behaviors.
    for (int color = 0; color < SimParameters::NUM_COLORS; color++) {
        for (int index = 0; index < this->boids[color].size(); index++) {
            Boid& us = this->boids[color][index];

            int num_friendly = 0;
            Eigen::Vector2d com = Vector2d::Zero(), vel = Vector2d::Zero(), avoid = Vector2d::Zero();

            // Examine all the friendly boids to gather relevant information.
            for (int ondex = 0; ondex < this->boids[color].size(); ondex++) {
                if (ondex == index) continue;

                Boid& them = this->boids[color][ondex];
                double distance = (them.position - us.position).norm();

                if (distance > params_->view_radius) continue;

                num_friendly++;
                com += them.position;
                vel += them.velocity;

                if (distance < params_->separation_threshold) avoid += (us.position - them.position);
            }

            // Then examine all the not friendly boids for other relevant information...
            for (int ocolor = 0; ocolor < SimParameters::NUM_COLORS; ocolor++) {
                if (ocolor == color) continue;
                
                for (int ondex = 0; ondex < this->boids[ocolor].size(); ondex++) {
                    Boid& them = this->boids[ocolor][ondex];
                    double distance = (them.position - us.position).norm();

                    if (distance > params_->view_radius) continue;
                    if (distance < params_->separation_threshold) avoid += (us.position - them.position);
                }
            }

            com /= num_friendly;
            vel /= num_friendly;

            Vector2d new_velocity = Vector2d::Zero();
            if (params_->separation_enabled) new_velocity += params_->separation_strength * avoid;
            if (params_->follow_com_enabled && num_friendly > 0) new_velocity += params_->follow_com_strength * (com - us.position);
            if (params_->follow_vel_enabled && num_friendly > 0) new_velocity += params_->follow_vel_strength * vel;
            if (params_->inertia_enabled) new_velocity += params_->inertia_strength * us.velocity;
            Vector2d norm_velocity = params_->speed * new_velocity.normalized();
            if (std::isnan(norm_velocity(0))) norm_velocity = Vector2d(params_->speed, 0.0);
            
            us.velocity = norm_velocity;
            us.position = us.position + us.velocity;
        }
    }

    return false;
}

// Spawn the given number of boids of the given color at the given location.
void BoidCore::spawn_boids(double x, double y, int count, SimParameters::Color color) {
    for (int i = 0; i < count; i++)
        this->boids[color].emplace_back(color, Vector2d(x, y), Vector2d::Zero());
}

// Create a boid spawner at the given location with the given color.
void BoidCore::create_spawner(double x, double y, SimParameters::Color color) {
    // TODO: implement me.
}

// Create a boid goal at the given location with the given color.
void BoidCore::create_goal(double x, double y, SimParameters::Color color) {
    // TODO: implement me.
}

std::tuple<MatrixXd, MatrixXi, MatrixXd> BoidCore::getCurrentMesh() const {
    // The shape of a boid; this points in the +X direction by default, and can be rotated.
    Vector2d BOID_SHAPE[3] = {
        Vector2d(0.03 - 0.03/3, 0), Vector2d(-0.03/3, 0.01), Vector2d(-0.03/3, -0.01)
    };

    // The RGB colors for each boid color.
    Vector4d COLOR_COLORS[4] = {
        Vector4d(1.0, 0.0, 0.0, 1.0),
        Vector4d(0.0, 1.0, 0.0, 1.0),
        Vector4d(1.0, 1.0, 0.0, 1.0),
        Vector4d(0.0, 0.0, 1.0, 1.0)
    };

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector4d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    // Render all of the boids based on the basic boid shape.
    for (int color = 0; color < SimParameters::NUM_COLORS; color++) {
        for (const Boid& boid : this->boids[color]) {
            double angle = std::atan2(boid.velocity(1), boid.velocity(0));
            if (std::isnan(angle)) angle = 0.0;
            Rotation2D<double> rot(angle);

            int start_index = verts.size();
            verts.push_back(z0(rot * BOID_SHAPE[0] + boid.position));
            vertexColors.push_back(COLOR_COLORS[color]);
            verts.push_back(z0(rot * BOID_SHAPE[1] + boid.position));
            vertexColors.push_back(COLOR_COLORS[color]);
            verts.push_back(z0(rot * BOID_SHAPE[2] + boid.position));
            vertexColors.push_back(COLOR_COLORS[color]);

            faces.emplace_back(start_index, start_index + 1, start_index + 2);
        }
    }

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(verts.size(), 3);
    renderC.resize(vertexColors.size(), 3);
    for (int i = 0; i < verts.size(); i++) {
        renderQ.row(i) = verts[i];
        renderC.row(i) = vertexColors[i];
    }
    renderF.resize(faces.size(), 3);
    for (int i = 0; i < faces.size(); i++)
        renderF.row(i) = faces[i];

    return std::make_tuple(renderQ, renderF, renderC);
}

}