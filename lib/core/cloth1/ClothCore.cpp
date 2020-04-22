#include "ClothCore.h"
#include "SimParameters.h"
#include <Eigen/SVD>
#include "igl/procrustes.h"
#include <iostream>
#include <limits>
#include <map>
#include <unordered_set>
#include <unordered_map>

using namespace Eigen;

namespace cloth1 {

/** A trivial hash implementation for a pair. */
struct hash_pair {
    template <class T, class U> size_t operator()(const std::pair<T, U>& p) const {
        return std::hash<T>{}(p.first) ^ std::hash<U>{}(p.second);
    }
};

ClothCore::ClothCore() {
    params_ = std::make_shared<SimParameters>();
    clickedVertex_ = -1;
}

ClothCore::~ClothCore() { }

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> ClothCore::getCurrentMesh() const {
    Eigen::MatrixXd C;
    return std::make_tuple(Q_, F_, C);
}

void ClothCore::attachMesh(Eigen::Ref<Eigen::MatrixX3d> V,
                           Eigen::Ref<Eigen::MatrixX3i> F,
                           double scale) {
    F_ = F;
    origQ_ = V * scale;
    initSimulation();
}

void ClothCore::initSimulation() {
    Q_ = origQ_;
    Qdot_.resize(Q_.rows(), 3);
    Qdot_.setZero();
    pinnedVerts_.clear();

    int nverts = Q_.rows();

    // Find the vertices closest to (1, 1, 0) and (-1, 1, 0) and treat them as our 'top left' and 'top right' vertices.
    int topleft = -1;
    int topright = -1;
    double topleftdist = -std::numeric_limits<double>::infinity();
    double toprightdist = -std::numeric_limits<double>::infinity();
    Eigen::Vector3d tr(1, 1, 0);
    Eigen::Vector3d tl(-1, 1, 0);
    for (int i = 0; i < nverts; i++) {
        double disttr = tr.dot(Q_.row(i));
        if (disttr > toprightdist) {
            toprightdist = disttr;
            topright = i;
        }
        double disttl = tl.dot(Q_.row(i));
        if (disttl > topleftdist) {
            topleftdist = disttl;
            topleft = i;
        }
    }
    pinnedVerts_.push_back(topleft);
    pinnedVerts_.push_back(topright);

    computeHinges();
}

bool ClothCore::simulateOneStep() {
    int nverts = Q_.rows();

    Eigen::MatrixXd oldQ = Q_;
    Q_ += params_->dt * Qdot_;

    // apply constraints
    for (int i = 0; i < params_->constraintIters; i++) applyConstraints();

    Qdot_ = (Q_ - oldQ) / params_->dt;

    if (params_->gravityEnabled)
        Qdot_.col(1) = Qdot_.col(1).array() + params_->dt * params_->gravityG;

    return false;
}

void ClothCore::hold(int vertex_id, const Eigen::Vector3d& position) {
    clickedVertex_ = vertex_id;
    curPos_ = position;
}

void ClothCore::updateHold(const Eigen::Vector3d& position) {
    if (clickedVertex_ < 0) return;
    curPos_ = position;
}

void ClothCore::releaseHold() {
    clickedVertex_ = -1;
}

void ClothCore::computeHinges() {
    // Iterate through all faces, and associate each face with it's edges (resulting in two faces/edge total for most edges).
    std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> edge_to_face;

    for (int face = 0; face < F_.rows(); face++) {
        int v1 = F_(face, 0), v2 = F_(face, 1), v3 = F_(face, 2);
        std::pair<int, int> edge1 = std::make_pair(std::min(v1, v2), std::max(v1, v2));
        std::pair<int, int> edge2 = std::make_pair(std::min(v1, v3), std::max(v1, v3));
        std::pair<int, int> edge3 = std::make_pair(std::min(v2, v3), std::max(v2, v3));

        edge_to_face[edge1].push_back(face);
        edge_to_face[edge2].push_back(face);
        edge_to_face[edge3].push_back(face);
    }

    // Count the total number of hinges so we can preallocate.
    int num_hinges = 0;
    for (auto& it : edge_to_face) {
        if (it.second.size() != 2) continue;
        num_hinges++;
    }

    H_.resize(num_hinges, 4);
    H_.setZero();
    int index = 0;

    // Then iterate through all the edges in the map, creating hinges out of ones that have two corresponding faces.
    for (auto& it : edge_to_face) {
        if (it.second.size() != 2) continue;

        std::unordered_set<int> vertices;
        vertices.insert(it.first.first);
        vertices.insert(it.first.second);

        int face1 = it.second[0], face2 = it.second[1];
        for (int index = 0; index < 3; index++) {
            vertices.insert(F_(face1, index));
            vertices.insert(F_(face2, index));
        }

        // We should end up with four distinct vertices; if not, give up.
        if (vertices.size() != 4) continue;

        auto iter = vertices.begin();
        H_(index, 0) = *(iter++);
        H_(index, 1) = *(iter++);
        H_(index, 2) = *(iter++);
        H_(index, 3) = *(iter++);
        index++;
    }
}

void ClothCore::applyConstraints() {
    // Interpolate pinned vertices towards rest locations.
    if (params_->pinEnabled) {
        double weight = params_->pinWeight;

        for (int index : pinnedVerts_)
            Q_.row(index) = weight * origQ_.row(index) + (1 - weight) * Q_.row(index);
    }

    // Rigidly map triangles in the rest configuration to triangles in the current configuration (this forces edge
    // lengths to remain the same in the mapping, limiting stretching).
    if (params_->stretchEnabled) {
        double weight = params_->stretchWeight;

        // For each triangle in the mesh, map it to the corresponding current mesh via the orthogonal Procrustes problem.
        for (int face = 0; face < F_.rows(); face++) {
            int v1 = F_(face, 0), v2 = F_(face, 1), v3 = F_(face, 2);

            Vector3d translation;
            Matrix3d orig_points, curr_points, rotation;
            orig_points << origQ_.row(v1), origQ_.row(v2), origQ_.row(v3);
            curr_points << Q_.row(v1), Q_.row(v2), Q_.row(v3);

            igl::procrustes(orig_points, curr_points, rotation, translation);

            Matrix3d transformed = (orig_points * rotation).rowwise() + translation.transpose();

            // Map the original points using the computed rotation and translation.
            Q_.row(v1) = weight * transformed.row(0) + (1 - weight) * curr_points.row(0);
            Q_.row(v2) = weight * transformed.row(1) + (1 - weight) * curr_points.row(1);
            Q_.row(v3) = weight * transformed.row(2) + (1 - weight) * curr_points.row(2);
        }
    }

    if (params_->bendingEnabled) {
        double weight = params_->bendingWeight;

        // For each hinge in the mesh, map it to the corresponding current mesh via the orthogonal Procrustes problem.
        for (int hinge = 0; hinge < H_.rows(); hinge++) {
            int v1 = H_(hinge, 0), v2 = H_(hinge, 1), v3 = H_(hinge, 2), v4 = H_(hinge, 3);

            Vector3d translation;
            Matrix3d rotation;
            Matrix<double, 4, 3> orig_points, curr_points;
            orig_points << origQ_.row(v1), origQ_.row(v2), origQ_.row(v3), origQ_.row(v4);
            curr_points << Q_.row(v1), Q_.row(v2), Q_.row(v3), Q_.row(v4);

            igl::procrustes(orig_points, curr_points, rotation, translation);

            Matrix<double, 4, 3> transformed = (orig_points * rotation).rowwise() + translation.transpose();

            // Map the original points using the computed rotation and translation.
            Q_.row(v1) = weight * transformed.row(0) + (1 - weight) * curr_points.row(0);
            Q_.row(v2) = weight * transformed.row(1) + (1 - weight) * curr_points.row(1);
            Q_.row(v3) = weight * transformed.row(2) + (1 - weight) * curr_points.row(2);
            Q_.row(v4) = weight * transformed.row(3) + (1 - weight) * curr_points.row(3);
        }
    }

    // Pin the clicked vertex at the mouse position.
    if (clickedVertex_ != -1 && params_->pullingEnabled) {
        double weight = params_->pullingWeight;

        Vector3d qPos = Q_.row(clickedVertex_);
        Q_.row(clickedVertex_) = weight * curPos_ + (1 - weight) * qPos;
    }
}

}
