#include "BoidCore.h"
#include "SimParameters.h"

using namespace Eigen;

namespace boid {

BoidCore::BoidCore() {
    params_ = std::make_shared<SimParameters>();
}

BoidCore::~BoidCore() { }

void BoidCore::initSimulation() {
    // TODO: implement me.
}

bool BoidCore::simulateOneStep() {
    // TODO: implement me.
    return false;
}

// Spawn the given number of boids of the given color at the given location.
void BoidCore::spawn_boids(double x, double y, int count, SimParameters::Color color) {
    // TODO: implement me.
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
    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

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