#include "BirdsCore.h"
#include "helper.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include <iostream>

#include <Eigen/LU> // Required for .inverse()
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

namespace bird1 {

/** A generic implementation of Newton's method. */
template<class F, class J> Vector3d newtons(Vector3d initial, int max_iters, double tolerance, const F&& func, const J&& jacobian) {
    Vector3d guess = initial;
    Vector3d value = func(initial);

    int steps = 0;
    while (steps < max_iters && value.norm() > tolerance) {
        auto jacob = jacobian(guess);

        guess = jacob.llt().solve(-value) + guess;
        value = func(guess);
        steps++;
    }

    std::cout << "Newton converged in " << steps << " steps (error " << value.norm() << ")" << std::endl;

    return guess;
}

BirdsCore::BirdsCore() {
    params_ = std::make_shared<SimParameters>();
    initSimulation();
}

BirdsCore::~BirdsCore() { }

std::tuple<MatrixXd, MatrixXi, MatrixXd> BirdsCore::getCurrentMesh() const {
    int totverts = 0;
    int totfaces = 0;
    for (const auto& rbi : bodies_) {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;
    for (const auto& rbi : bodies_) {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++) {
            for (int j = 0; j < 3; j++) {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }

    return std::make_tuple(renderQ, renderF, renderC);
}


void BirdsCore::initSimulation() {
    rigid_body_id_ = 0;
    time_ = 0;
}

bool BirdsCore::simulateOneStep() {
    time_ += params_->timeStep;
    int nbodies = (int)bodies_.size();
    double h = params_->timeStep;

    // The total force on every center of mass; the segment (3*i, 3) corresponds to the force on the ith body.
    VectorXd com_force = this->force_c();

    // A velocity-verlet looking time integrator which updates four variables: next center of masses, next orientations,
    // next velocities, and next angular velocities. Since the bodies do not interact in any way, we simulate each one
    // separately for now. This will change in the next milestone (though the code should be the same, just done on
    // larger vectors).
    for (int i = 0; i < nbodies; i++) {
        RigidBodyInstance& body = *bodies_[i];

        // Scale the inertia tensor by this specific instance density.
        Matrix3d inertia_tensor = body.getTemplate().getInertiaTensor() * body.density;
        // We only actually need the inverse mass matrix, which we will compute explicitly as 1/(density*volume).
        double mass = body.getTemplate().getVolume() * body.density;
        Matrix3d inv_mass;
        inv_mass <<
            1./mass, 0, 0,
            0, 1./mass, 0,
            0, 0, 1./mass;

        // Easy part: step the center of mass and rotation by their velocities directly.
        Vector3d next_c = body.c + h * body.cvel;
        Vector3d next_theta = VectorMath::axisAngle(VectorMath::rotationMatrix(body.theta) * VectorMath::rotationMatrix(h * body.w));

        // Extract the force acting on this specific body. There is currently no rotational force anywhere.
        Vector3d c_force = com_force.segment(3*i, 3);
        Vector3d theta_force = Vector3d::Zero();
        
        // Medium part: compute the next velocity explicitly.
        Vector3d next_cvel = body.cvel - h * inv_mass * c_force;

        // Hard part: compute the next angular velocity implicitly using newton's method.
        Vector3d next_w = newtons(body.w, params_->NewtonMaxIters, params_->NewtonTolerance,
            // The function itself.
            [&](const Vector3d& w) {
                return w.transpose() * inertia_tensor * VectorMath::TMatrix(-h * w).inverse()
                    - body.w.transpose() * inertia_tensor * VectorMath::TMatrix(h * body.w).inverse()
                    + h * theta_force.transpose() * VectorMath::TMatrix(next_theta).inverse();
            },
            // The (quasi-)derivative of the function; the full derivative is annoying, so we only implement the
            // dominant term (the other term is a factor of h smaller which is negligible for h ~ 10^-6)
            [&](const Vector3d& w) {
                // TODO: I am very braindead, the t-matrix may need to be transposed (the inertia tensor is symmetric
                // and does not need to be).
                return VectorMath::TMatrix(-h * w).inverse().transpose() * inertia_tensor;
                // return inertia_tensor * VectorMath::TMatrix(-h * w).inverse().transpose();
            });

        // And persist our computed values.
        body.c = next_c;
        body.theta = next_theta;
        body.cvel = next_cvel;
        body.w = next_w;
    }

    return false;
}

VectorXd BirdsCore::force_c() {
    if (params_->gravityEnabled) return this->force_gravity();
    else return VectorXd::Zero(3*bodies_.size());
}

VectorXd BirdsCore::force_gravity() {
    VectorXd result = VectorXd::Zero(3*bodies_.size());

    for (int i = 0; i < bodies_.size(); i++) {
        for (int j = i + 1; j < bodies_.size(); j++) {
            Vector3d c1 = bodies_[i]->c;
            Vector3d c2 = bodies_[j]->c;
            double dist = (c1 - c2).norm();

            Vector3d force = params_->gravityG * (bodies_[i]->mass() * bodies_[j]->mass()) / (dist * dist) * (c1 - c2).normalized();

            result.segment(3*i, 3) += force;
            result.segment(3*j, 3) -= force;
        }
    }

    return result;
}

void BirdsCore::clearScene() {
    bodies_.clear();
    templates_.clear();
    init_configurations_.clear();
}

VectorXi BirdsCore::addMesh(const std::string& file_name, double scale, const Eigen::MatrixXd& Qs) {
    auto tup = loadOBJ(file_name);
    templates_.emplace_back(new RigidBodyTemplate(std::get<0>(tup), std::get<1>(tup), scale));

    auto rbt = templates_.back();
    init_configurations_.emplace_back(Qs);
    return addInstances(rbt, Qs);
}

std::shared_ptr<RigidBodyInstance> BirdsCore::queryRigidBodyInstance(int32_t bid) {
    for (auto& b : bodies_)
        if (b->bid == bid) return b;
    return std::shared_ptr<RigidBodyInstance>(nullptr);
}

int32_t BirdsCore::addSingleInstance(std::shared_ptr<RigidBodyTemplate> rbt,
                             double density,
                             const Eigen::Vector3d &c,
                             const Eigen::Vector3d &theta,
                             const Eigen::Vector3d &cvel,
                             const Eigen::Vector3d &w) {
    bodies_.emplace_back(new RigidBodyInstance(*rbt, c, theta, cvel, w, density));
    bodies_.back()->bid = rigid_body_id_++;
    return bodies_.back()->bid;
}

VectorXi BirdsCore::addInstances(std::shared_ptr<RigidBodyTemplate> rbt, const Eigen::MatrixXd& Qs) {
    Eigen::VectorXi ret;
    ret.resize(Qs.rows());
    for (int i = 0; i < Qs.rows(); i++) {
        double density;
        Eigen::Vector3d c, cvel, theta, w;
        density = Qs(i, 0);
        int base = 1;
        c << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        cvel << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        theta << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        w << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        ret(i) = addSingleInstance(rbt, density, c, theta, cvel, w);
    }
    return ret;
}

}
