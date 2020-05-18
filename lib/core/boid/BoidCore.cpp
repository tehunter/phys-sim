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

template<class T> constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

BoidCore::BoidCore() {
    params_ = std::make_shared<SimParameters>();
}

BoidCore::~BoidCore() { }

void BoidCore::initSimulation() {
    this->boids.clear();
    this->predators.clear();
    this->obstacles.clear();
    this->grid.clear();

    // Initialize the grid layout to have (NUM_GRID_CELLS x NUM_GRID_CELLS).
    for (int y = 0; y < SimParameters::NUM_GRID_CELLS; y++) {
        this->grid.emplace_back();
        this->grid[y].resize(SimParameters::NUM_GRID_CELLS);
    }

    for (int color = 0; color < SimParameters::NUM_COLORS; color++)
        this->goals[color] = std::make_pair(false, Vector2d::Zero());

    // Place down some initial boids.
    this->spawn_boids(0.0, 0.0, 1, SimParameters::RED);
    this->spawn_boids(0.05, 0.0, 1, SimParameters::RED);
    this->spawn_boids(-0.05, 0.0, 1, SimParameters::RED);
    this->spawn_boids(0.0, 0.05, 1, SimParameters::RED);
    this->spawn_boids(0.0, -0.05, 1, SimParameters::RED);
}

bool BoidCore::simulateOneStep() {
    // Clear the grid so that we can re-place boids in it.
    for (int r = 0; r < SimParameters::NUM_GRID_CELLS; r++) {
        for (int c = 0; c < SimParameters::NUM_GRID_CELLS; c++)
            this->grid[r][c].clear();
    }

    // Index the boids into the grid for fast spatial lookups.
    for (int index = 0; index < this->boids.size(); index++) {
        Boid& us = this->boids[index];
        Vector2i loc = this->grid_cell(us.position);
        this->grid[loc(0)][loc(1)].push_back(index);
    }

    // Update all boid velocities...
    std::vector<int> neighbors;
    for (int index = 0; index < this->boids.size(); index++) {
        Boid& us = this->boids[index];

        if (this->fill_neighbors(us.position, params_->view_radius, index, neighbors)) {
            if (params_->follow_com_enabled) us.velocity += params_->follow_com_strength * params_->speed * this->velocity_com(us, neighbors).normalized();
            if (params_->follow_vel_enabled) us.velocity += params_->follow_vel_strength * this->velocity_vel(us, neighbors);
            if (params_->separation_enabled) us.velocity += params_->separation_strength * params_->speed * this->velocity_avoid(us, neighbors).normalized();

            neighbors.clear();
        }

        // Handle the edges of the arena, goals, and predators.
        if (params_->walls_enabled) us.velocity += params_->wall_strength * params_->speed * this->velocity_wall(us.position);
        if (params_->goals_enabled) us.velocity += params_->goal_strength * params_->speed * this->velocity_goal(us).normalized();
        if (params_->avoid_predators_enabled) us.velocity += params_->avoid_predators_strength * params_->speed * this->velocity_avoid_preds(us).normalized();

        us.velocity = us.velocity.normalized() * params_->speed;
        us.position += us.velocity;
    }

    // Then update predator velocities...
    for (int index = 0; index < this->predators.size(); index++) {
        Predator& us = this->predators[index];

        if (this->fill_neighbors(us.position, params_->predator_view_radius * (us.size_factor() - Predator::BASE_SIZE_FACTOR + 1), -1, neighbors)) {
            // Find the closest neighbor boid, and aggressively try to eat it.
            int closest = -1;
            double distance = 99999999;
            for (int ondex : neighbors) {
                double d = (this->boids[ondex].position - us.position).norm();
                if (d >= distance) continue;

                closest = ondex;
                distance = d;
            }

            us.velocity += params_->predator_chase_strength * params_->predator_speed * (this->boids[closest].position - us.position).normalized();
            
            // If the boid is right next to us, make it die.
            if (distance < us.size_factor() * 0.025) {
                this->boids[closest].dead = true;
                us.eaten += 1;
            }

            neighbors.clear();
        }

        // Handle the edges of the arena, and goals.
        if (params_->walls_enabled) us.velocity += params_->wall_strength * params_->predator_speed * this->velocity_wall(us.position);

        us.velocity = us.velocity.normalized() * params_->predator_speed;
        us.position += us.velocity;
    }

    // Clean up dead boids.
    for (int index = this->boids.size() - 1; index >= 0; index--) {
        if (this->boids[index].dead) this->boids.erase(this->boids.begin() + index);
    }

    return false;
}

// Spawn the given number of boids of the given color at the given location.
void BoidCore::spawn_boids(double x, double y, int count, SimParameters::Color color) {
    for (int i = 0; i < count; i++)
        this->boids.emplace_back(color, Vector2d(x, y) + this->random_initial_velocity() * 20, random_initial_velocity());
}

void BoidCore::spawn_predator(double x, double y) {
    this->predators.emplace_back(Vector2d(x, y), random_initial_velocity());
}

void BoidCore::create_obstacle(double x, double y) {
    // TODO: Make radius configurable.
    this->obstacles.emplace_back(Vector2d(x, y), 0.05);
}

// Create a boid goal at the given location with the given color.
void BoidCore::create_goal(double x, double y, SimParameters::Color color) {
    // If a goal exists and we click close to it, then remove it.
    Vector2d click = Vector2d(x, y);

    if (this->goals[color].first && (this->goals[color].second - click).norm() < SimParameters::GOAL_REMOVE_DIST) {
        this->goals[color] = std::make_pair(false, Vector2d::Zero());
    } else {
        this->goals[color] = std::make_pair(true, click);
    }
}

void BoidCore::update_mouse_position(double x, double y) {
    this->mouse_location = Vector2d(x, y);
}

Vector2d BoidCore::random_initial_velocity() {
    std::uniform_real_distribution<double> dist(-M_PI, M_PI);
    double angle = dist(this->rng);
    return params_->speed * Vector2d(std::cos(angle), std::sin(angle));
}

Vector2i BoidCore::grid_cell(Vector2d position) const {
    double x = clamp(position(0), SimParameters::GRID_START_X, SimParameters::GRID_END_X);
    double y = clamp(position(1), SimParameters::GRID_START_Y, SimParameters::GRID_END_Y);

    int row = clamp((int) ((y - SimParameters::GRID_START_Y) / SimParameters::GRID_CELL_HEIGHT), 0, SimParameters::NUM_GRID_CELLS - 1);
    int col = clamp((int) ((x - SimParameters::GRID_START_X) / SimParameters::GRID_CELL_WIDTH), 0, SimParameters::NUM_GRID_CELLS - 1);
    
    return Vector2i(row, col);
}

int BoidCore::fill_neighbors(Vector2d position, double radius, int ignore, std::vector<int>& output) const {
    Vector2i grid_pos = grid_cell(position);

    int xSearch = (int) std::ceil(radius / params_->GRID_CELL_WIDTH);
    int ySearch = (int) std::ceil(radius / params_->GRID_CELL_HEIGHT);

    int count = 0;
    for (int y = -ySearch + grid_pos[0]; y <= ySearch + grid_pos[0]; y++) {
        for (int x = -xSearch + grid_pos[1]; x <= xSearch + grid_pos[1]; x++) {
            if (x < 0 || x >= grid[y].size() || y < 0 || y >= grid.size()) continue;

            for (int index : this->grid[y][x]) {
                if (index == ignore) continue;

                double distance = (this->boids[index].position - position).norm();
                if (distance > radius) continue;

                output.push_back(index);
                count++;
            }
        }
    }

    return count;
}

std::pair<bool, Vector2d> BoidCore::goal(SimParameters::Color color) const {
    // If there is a goal, prefer that over everything else.
    if (this->goals[color].first) return this->goals[color];

    // If the mouse mode is 'follow', then prefer to follow around the mouse.
    if (params_->mouse_mode == SimParameters::MM_FOLLOW) return std::make_pair(true, this->mouse_location);

    // Otherwise, no goal.
    return std::make_pair(false, Vector2d::Zero());
}

Vector2d BoidCore::velocity_com(const Boid& us, std::vector<int>& neighbors) const {
    int count = 0;
    Vector2d com = Vector2d::Zero();

    for (int neighbor : neighbors) {
        if (this->boids[neighbor].color != us.color) continue;

        com += this->boids[neighbor].position;
        count++;
    }

    if (count == 0) return Vector2d::Zero();
    return (com / (float) count) - us.position;
}

Vector2d BoidCore::velocity_vel(const Boid& us, std::vector<int>& neighbors) const {
    int count = 0;
    Vector2d velocity = Vector2d::Zero();

    for (int neighbor : neighbors) {
        if (this->boids[neighbor].color != us.color) continue;

        velocity += this->boids[neighbor].velocity;
        count++;
    }

    if (count == 0) return Vector2d::Zero();
    return (velocity / (float) count) - us.velocity;
}

Vector2d BoidCore::velocity_avoid(const Boid& us, std::vector<int>& neighbors) const {
    Vector2d avoid = Vector2d::Zero();

    for (int neighbor : neighbors) {
        Vector2d diff = (us.position - this->boids[neighbor].position);
        if (diff.norm() > params_->separation_threshold) continue;

        avoid += diff;
    }

    return avoid;
}

Vector2d BoidCore::velocity_wall(Vector2d position) const {
    Vector2d result(0.0, 0.0);
    if (position(0) < -SimParameters::WALL_X + SimParameters::WALL_MARGIN) result(0) += 1.0;
    else if (position(0) > SimParameters::WALL_X - SimParameters::WALL_MARGIN) result(0) -= 1.0;

    if (position(1) < -SimParameters::WALL_Y + SimParameters::WALL_MARGIN) result(1) += 1.0;
    else if (position(1) > SimParameters::WALL_Y - SimParameters::WALL_MARGIN) result(1) -= 1.0;

    return result.normalized();
}

Vector2d BoidCore::velocity_goal(const Boid& us) const {
    // Check if we have a goal in the first place.
    auto goal = this->goal(us.color);
    if (!goal.first) return Vector2d::Zero();

    return (goal.second - us.position);
}

Vector2d BoidCore::velocity_avoid_preds(const Boid& us) const {
    Vector2d avoid = Vector2d::Zero();

    for (const Predator& neighbor : this->predators) {
        if ((neighbor.position - us.position).norm() > params_->view_radius + 0.025 * neighbor.size_factor()) continue;
        avoid += (us.position - neighbor.position);
    }

    return avoid;
}

std::tuple<MatrixXd, MatrixXi, MatrixXd> BoidCore::getCurrentMesh() const {
    // The shape of a boid; this points in the +X direction by default, and can be rotated.
    Vector2d BOID_SHAPE[3] = {
        Vector2d(0.02 - 0.02/3, 0), Vector2d(-0.02/3, 0.007), Vector2d(-0.02/3, -0.007)
    };

    // The RGB colors for each boid color.
    Vector4d COLOR_COLORS[4] = {
        Vector4d(1.0, 0.0, 0.0, 1.0),
        Vector4d(0.0, 1.0, 0.0, 1.0),
        Vector4d(1.0, 1.0, 0.0, 1.0),
        Vector4d(0.0, 0.0, 1.0, 1.0)
    };

    // The color of a predator (white).
    Vector4d PREDATOR_COLOR = Vector4d(1.0, 1.0, 1.0, 1.0);
    // The color of an obstacle (black).
    Vector4d OBSTACLE_COLOR = Vector4d(0.0, 0.0, 0.0, 1.0);

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector4d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    // Render all of the boids based on the basic boid shape.
    for (const Boid& boid : this->boids) {
        double angle = std::atan2(boid.velocity(1), boid.velocity(0));
        if (std::isnan(angle)) angle = 0.0;
        Rotation2D<double> rot(angle);

        int start_index = verts.size();
        verts.push_back(z0(rot * BOID_SHAPE[0] + boid.position));
        vertexColors.push_back(COLOR_COLORS[boid.color]);
        verts.push_back(z0(rot * BOID_SHAPE[1] + boid.position));
        vertexColors.push_back(COLOR_COLORS[boid.color]);
        verts.push_back(z0(rot * BOID_SHAPE[2] + boid.position));
        vertexColors.push_back(COLOR_COLORS[boid.color]);

        faces.emplace_back(start_index, start_index + 1, start_index + 2);
    }

    // Render predators, which can get very chonky and big.
    for (const Predator& predator : this->predators) {
        double angle = std::atan2(predator.velocity(1), predator.velocity(0));
        if (std::isnan(angle)) angle = 0.0;
        Rotation2D<double> rot(angle);

        int start_index = verts.size();
        verts.push_back(z0(rot * BOID_SHAPE[0] * predator.size_factor() + predator.position));
        vertexColors.push_back(PREDATOR_COLOR);
        verts.push_back(z0(rot * BOID_SHAPE[1] * predator.size_factor() + predator.position));
        vertexColors.push_back(PREDATOR_COLOR);
        verts.push_back(z0(rot * BOID_SHAPE[2] * predator.size_factor() + predator.position));
        vertexColors.push_back(PREDATOR_COLOR);

        faces.emplace_back(start_index, start_index + 1, start_index + 2);
    }

    // Render goals with the appropriate color.
    for (int color = 0; color < SimParameters::NUM_COLORS; color++) {
        auto goal = this->goals[color];
        if (!goal.first) continue;

        int start_index = verts.size();
        verts.push_back(z0(goal.second + Vector2d(0, 0.03)));
        vertexColors.push_back(COLOR_COLORS[color]);
        verts.push_back(z0(goal.second + Vector2d(-0.03, 0)));
        vertexColors.push_back(COLOR_COLORS[color]);
        verts.push_back(z0(goal.second + Vector2d(0, -0.03)));
        vertexColors.push_back(COLOR_COLORS[color]);
        verts.push_back(z0(goal.second + Vector2d(0.03, 0)));
        vertexColors.push_back(COLOR_COLORS[color]);

        faces.emplace_back(start_index, start_index + 1, start_index + 2);
        faces.emplace_back(start_index, start_index + 2, start_index + 3);
    }

    // Render obstacles as circle-looking things.

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(verts.size(), 3);
    renderC.resize(vertexColors.size(), 4);
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