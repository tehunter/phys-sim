#include "BirdsCore.h"
#include "helper.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "CollisionDetection.h"
#include <unordered_map>
#include <iostream>
#include <Eigen/LU> // Required for .inverse()
#include <Eigen/Geometry> // Required for .cross()

using namespace Eigen;

namespace bird2 {

/** A trivial hash implementation for a pair. */
struct hash_pair {
    template <class T, class U> size_t operator()(const std::pair<T, U>& p) const {
        return std::hash<T>{}(p.first) ^ std::hash<U>{}(p.second);
    }
};

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

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
BirdsCore::getCurrentMesh() const {
    int totverts = 0;
    int totfaces = 0;

    // floor
    totverts += 5;
    totfaces += 4;

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

    // floor
    double floory = -1.0;
    renderQ.row(0) << 0, floory, 0;
    renderQ.row(1) << 1e6, floory, 1e6;
    renderQ.row(2) << -1e6, floory, 1e6;
    renderQ.row(3) << -1e6, floory, -1e6;
    renderQ.row(4) << 1e6, floory, -1e6;
    voffset += 5;

    renderF.row(0) << 0, 2, 1;
    renderF.row(1) << 0, 3, 2;
    renderF.row(2) << 0, 4, 3;
    renderF.row(3) << 0, 1, 4;
    foffset += 4;

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
    double h = params_->timeStep;
    int nbodies = (int)bodies_.size();

    // Step forward center of mass/rotation according to velocities, since velocity verlet computes forces using the
    // next timesteps COM/theta.
    std::vector<Vector3d> old_theta(bodies_.size());
    for (int i = 0; i < nbodies; i++) {
        RigidBodyInstance& body = *bodies_[i];

        // Easy part: step the center of mass and rotation by their velocities directly.
        Vector3d next_c = body.c + h * body.cvel;
        Vector3d next_theta = VectorMath::axisAngle(VectorMath::rotationMatrix(body.theta) * VectorMath::rotationMatrix(h * body.w));
        
        old_theta[i] = body.theta;

        body.c = next_c;
        body.theta = next_theta;
    }

    time_ += h;
    std::set<Collision> collisions = collisionDetection(bodies_);

    VectorXd com_force = VectorXd::Zero(3 * nbodies);
    VectorXd theta_force = VectorXd::Zero(3 * nbodies);

    // Compute the forces acting on all bodies (both translational AND rotational).
    this->computeForces(com_force, theta_force);
    this->computePenaltyCollisionForces(collisions, com_force, theta_force);

    // Apply impulse forces due to collisions directly. These modify body velocities instantaneously, and so
    // are not part of the time integration step below.
    this->applyCollisionImpulses(collisions);

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

        // Extract the force acting on this specific body. There is currently no rotational force anywhere.
        Vector3d local_com_force = com_force.segment(3*i, 3);
        Vector3d local_theta_force = theta_force.segment(3*i, 3);
        
        // Medium part: compute the next velocity explicitly.
        Vector3d next_cvel = body.cvel + h * inv_mass * local_com_force;

        // Hard part: compute the next angular velocity implicitly using newton's method.
        Vector3d next_w = newtons(body.w, params_->NewtonMaxIters, params_->NewtonTolerance,
            // The function itself.
            [&](const Vector3d& w) {
                // TODO: Use oldtheta instead of current body theta.
                return w.transpose() * inertia_tensor * VectorMath::TMatrix(-h * w).inverse()
                    - body.w.transpose() * inertia_tensor * VectorMath::TMatrix(h * body.w).inverse()
                    - h * local_theta_force.transpose() * VectorMath::TMatrix(old_theta[i]).inverse();
            },
            // The (quasi-)derivative of the function; the full derivative is annoying, so we only implement the
            // dominant term (the other term is a factor of h smaller which is negligible for h ~ 10^-6)
            [&](const Vector3d& w) {
                // TODO: May need to be transposed?
                return VectorMath::TMatrix(-h * w).inverse().transpose() * inertia_tensor;
            });

        // And persist our computed values.
        body.cvel = next_cvel;
        body.w = next_w;
    }

    return false;
}

void BirdsCore::clearScene() {
    bodies_.clear();
    templates_.clear();
    init_configurations_.clear();
}

Eigen::VectorXi BirdsCore::addMesh(const std::string& file_name,
                   double scale, const Eigen::MatrixXd& Qs) {
    auto tup = bird1::loadOBJ(file_name);
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

Eigen::VectorXi BirdsCore::addInstances(std::shared_ptr<RigidBodyTemplate> rbt,
                        const Eigen::MatrixXd& Qs) {
    Eigen::VectorXi ret;
    ret.resize(Qs.rows());
    for (int i = 0; i < Qs.rows(); i++) {
        double density;
        Eigen::Vector3d c, cvel, theta, w;
        density = Qs(i, 0);
        int base = 1;
        c << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        theta << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        cvel << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        w << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        ret(i) = addSingleInstance(rbt, density, c, theta, cvel, w);
    }
    return ret;
}

void BirdsCore::computeForces(VectorXd &Fc, VectorXd &Ftheta) {
    if (!params_->gravityEnabled) return;

    // Apply gravitational force to the COM of all bodies.
    for (int i = 0; i < bodies_.size(); i++) {
        double mass = bodies_[i]->getTemplate().getVolume() * bodies_[i]->density;
        Fc[3 * i + 1] -= params_->gravityG * mass;
    }
}

void BirdsCore::computePenaltyCollisionForces(const std::set<Collision>& collisions,
                                              Eigen::Ref<VectorXd> Fc,
                                              Eigen::Ref<VectorXd> Ftheta) {
    if (!params_->penaltyEnabled) return;

    for (const Collision& collision : collisions) {
        // Here, a body1 vertex is colliding with a tetrahedron in body2 (or the floor, if collision.body2 == -1).
        // Extract the first body center of mass & orientation, as well as the colliding vertex IN BODY1 TEMPLATE
        // COORDINATES.
        Vector3d com1 = bodies_[collision.body1]->c;
        Vector3d theta1 = bodies_[collision.body1]->theta;
        Vector3d vertex1 = bodies_[collision.body1]->getTemplate().getVerts().row(collision.collidingVertex);

        // We need to specially handle the case when the colliding object is the floor, since the distance function used
        // is different. Otherwise, we use the normal body distance function.
        if (collision.body2 == -1) {
            // The vertex in 'floor coordinates' (where y=0 is the floor).
            Vector3d vertex2 = VectorMath::rotationMatrix(theta1) * vertex1 + com1 - Vector3d(0, -1, 0);

            // The distance function is just the y coordinate of the vertex.
            Vector3d shared = -params_->penaltyStiffness * vertex2(1) * Vector3d(0, 1, 0);

            // Now the black magic math, derived from the potential function in the handout.
            Fc.segment<3>(3*collision.body1) += shared.transpose();
            Ftheta.segment<3>(3*collision.body1) -= shared.transpose() * VectorMath::rotationMatrix(theta1)
                * VectorMath::crossProductMatrix(vertex1) * VectorMath::TMatrix(theta1);
        } else {
            // Get the tetrahedra, center of mass, and orientation of the second body.
            Vector3d com2 = bodies_[collision.body2]->c, theta2 = bodies_[collision.body2]->theta;
            const RigidBodyTemplate& template2 = bodies_[collision.body2]->getTemplate();

            // The vertex in the body 2 template coordinates.
            Vector3d vertex2 = VectorMath::rotationMatrix(-theta2) * (VectorMath::rotationMatrix(theta1) * vertex1 + com1 - com2);
            
            // The shared multiplier for all the the forces.
            // TODO: Check the sign here, may not work.
            Vector3d shared = -params_->penaltyStiffness * template2.distance(vertex2, collision.collidingTet)
                * template2.Ddistance(collision.collidingTet);

            // Now the black magic math, derived from the potential function in the handout.
            Fc.segment<3>(3*collision.body1) += shared.transpose() * VectorMath::rotationMatrix(-theta2);
            Fc.segment<3>(3*collision.body2) -= shared.transpose() * VectorMath::rotationMatrix(-theta2);
            Ftheta.segment<3>(3*collision.body1) -= shared.transpose() * VectorMath::rotationMatrix(-theta2) * VectorMath::rotationMatrix(theta1)
                * VectorMath::crossProductMatrix(vertex1) * VectorMath::TMatrix(theta1);
            Ftheta.segment<3>(3*collision.body2) += shared.transpose() * VectorMath::rotationMatrix(-theta2)
                * VectorMath::crossProductMatrix(VectorMath::rotationMatrix(theta1) * vertex1 + com1 - com2)
                * VectorMath::TMatrix(-theta2);
        }
    }
}

/** Translate the point in the 'from' rigid body instance coordinate space to the one in the 'towards' coordinate space.*/
Vector3d translate_to(const RigidBodyInstance& from, Vector3d point, const RigidBodyInstance& towards) {
    return VectorMath::rotationMatrix(-towards.theta) * (VectorMath::rotationMatrix(from.theta) * point + from.c - towards.c);
}

/** Compute the relative velocity of the point in instance 'from' the point in instance 'towards'. */
double relative_velocity(const RigidBodyInstance& from, Vector3d point, const RigidBodyInstance& towards, int tet) {
    Vector3d inside = VectorMath::rotationMatrix(from.theta) * point + from.c - towards.c;
    Vector3d d_tr_point = -VectorMath::rotationMatrix(-towards.theta) * VectorMath::crossProductMatrix(towards.w) * inside
        + VectorMath::rotationMatrix(-towards.theta) * (VectorMath::rotationMatrix(from.theta) * VectorMath::crossProductMatrix(from.w) * point
            + from.cvel - towards.cvel);

    // The result is a 1x1 vector, extract value via sum.
    return towards.getTemplate().Ddistance(tet).dot(d_tr_point);
}

double floor_velocity(const RigidBodyInstance& from, Vector3d point) {
    return Vector3d(0, -1, 0).dot(VectorMath::rotationMatrix(from.theta) * VectorMath::crossProductMatrix(from.w) * point + from.cvel);
}

void BirdsCore::applyCollisionImpulses(const std::set<Collision>& collisions) {
    if (!params_->impulsesEnabled) return;

    // For each pair of colliding bodies, find the (vertex, tet) pair with negative relative velocity and maximal
    // negative distance.
    std::unordered_map<std::pair<int, int>, double, hash_pair> distances;
    std::unordered_map<std::pair<int, int>, Collision, hash_pair> best_collisions;
    for (const auto& collision : collisions) {
        Vector3d vertex1 = bodies_[collision.body1]->getTemplate().getVerts().row(collision.collidingVertex);

        double distance, velocity;
        // Specially handle collisions with the floor.
        if (collision.body2 == -1) {
            const RigidBodyInstance& body = *bodies_[collision.body1];
            Vector3d vertex2 = VectorMath::rotationMatrix(body.theta) * vertex1 + body.c - Vector3d(0, -1, 0);

            velocity = floor_velocity(*bodies_[collision.body1], vertex1);
            distance = -vertex2(1);
        } else {
            Vector3d vertex2 = translate_to(*bodies_[collision.body1], vertex1, *bodies_[collision.body2]);

            velocity = relative_velocity(*bodies_[collision.body1], vertex1, *bodies_[collision.body2], collision.collidingTet);
            distance = bodies_[collision.body2]->getTemplate().distance(vertex2, collision.collidingTet);
        }

        // Object is moving away, 
        if (velocity < 0.0 || distance < 0.0) continue;

        std::pair<int, int> key = std::make_pair(std::min(collision.body1, collision.body2), std::max(collision.body1, collision.body2));

        if (distances.count(key) == 0 || distances[key] < distance) {
            distances[key] = distance;
            best_collisions[key] = collision;
        }
    }

    for (auto const& cpair : best_collisions) {
        const Collision& coll = cpair.second;

        Vector3d vertex1 = bodies_[coll.body1]->getTemplate().getVerts().row(coll.collidingVertex);

        // Once again, specially handle floor vs. normal collisions.
        if (coll.body2 == -1) {
            RigidBodyInstance& body = *bodies_[coll.body1];

            // To find the direction to push the object back out of, we compute -[dg] (i.e., the derivative of the
            // constraint function).
            Vector3d c_dir(0, 1, 0);
            Vector3d theta_dir = Vector3d(0, -1, 0).transpose() * VectorMath::rotationMatrix(body.theta)
                * VectorMath::crossProductMatrix(vertex1) * VectorMath::TMatrix(body.theta);

            double velocity = floor_velocity(*bodies_[coll.body1], vertex1);
            
            // Scale the inertia tensor by this specific instance density.
            Matrix3d inertia_tensor = body.getTemplate().getInertiaTensor() * body.density;
            double mass = body.getTemplate().getVolume() * body.density;
            Matrix3d inv_mass_matrix = Matrix3d::Identity() * (1./mass);

            // Now, compute the alpha which reflects the appropriate impulse magnitude...
            double z = Vector3d(0, -1, 0).transpose() *
                (VectorMath::rotationMatrix(body.theta) * (VectorMath::crossProductMatrix(inertia_tensor.inverse().transpose()
                    * VectorMath::TMatrix(body.theta).inverse().transpose() * theta_dir) * vertex1)
                + inv_mass_matrix * c_dir);
            double alpha = -(1 + params_->CoR) * velocity / z;

            if (z == 0.0) continue;

            body.cvel += alpha * inv_mass_matrix * c_dir;
            body.w += alpha * inertia_tensor.inverse().transpose() * VectorMath::TMatrix(body.theta).inverse().transpose() * theta_dir;
        } else {
            RigidBodyInstance& body1 = *bodies_[coll.body1];
            RigidBodyInstance& body2 = *bodies_[coll.body2];

            Vector3d ddist = body1.getTemplate().Ddistance(coll.collidingTet);

            // To find the direction to push the object back out of, we compute -[dg] (i.e., the derivative of the
            // constraint function).
            Vector3d c1_dir = ddist.transpose() * VectorMath::rotationMatrix(-body2.theta);
            Vector3d c2_dir = -ddist.transpose() * VectorMath::rotationMatrix(-body2.theta);
            Vector3d theta1_dir = -ddist.transpose() * VectorMath::rotationMatrix(-body2.theta) * VectorMath::rotationMatrix(body1.theta)
                * VectorMath::crossProductMatrix(vertex1) * VectorMath::TMatrix(body1.theta);
            Vector3d theta2_dir = ddist.transpose() * VectorMath::rotationMatrix(-body2.theta) * VectorMath::crossProductMatrix(
                VectorMath::rotationMatrix(body1.theta) * vertex1 + body1.c - body2.c) * VectorMath::TMatrix(-body2.theta);

            // Some more constants we need.
            Matrix3d inertia_tensor1 = body1.getTemplate().getInertiaTensor() * body1.density;
            Matrix3d inertia_tensor2 = body2.getTemplate().getInertiaTensor() * body2.density;
            double mass1 = body1.getTemplate().getVolume() * body1.density;
            Matrix3d inv_mass_matrix1 = Matrix3d::Identity() * (1./mass1);
            double mass2 = body2.getTemplate().getVolume() * body2.density;
            Matrix3d inv_mass_matrix2 = Matrix3d::Identity() * (1./mass2);

            double velocity = relative_velocity(*bodies_[coll.body1], vertex1, *bodies_[coll.body2], coll.collidingTet);

            double z = ddist.transpose() * VectorMath::rotationMatrix(-body2.theta) * (VectorMath::rotationMatrix(body1.theta)
                * VectorMath::crossProductMatrix(inertia_tensor1.inverse().transpose() * VectorMath::TMatrix(body1.theta).inverse().transpose() * theta1_dir)
                * vertex1
                + inv_mass_matrix1 * c1_dir - inv_mass_matrix2 * c2_dir
                - VectorMath::crossProductMatrix(inertia_tensor2.inverse().transpose() * VectorMath::TMatrix(body2.theta).inverse().transpose() * theta2_dir)
                * (VectorMath::rotationMatrix(body1.theta) * vertex1 + body1.c - body2.c));
            double alpha = -(1 + params_->CoR) * velocity / z;

            if (z == 0.0) continue;

            body1.cvel += alpha * inv_mass_matrix1 * c1_dir;
            body2.cvel += alpha * inv_mass_matrix2 * c2_dir;
            body1.w += alpha * inertia_tensor1.inverse().transpose() * VectorMath::TMatrix(body1.theta).inverse().transpose() * theta1_dir;
            body2.w += alpha * inertia_tensor2.inverse().transpose() * VectorMath::TMatrix(body2.theta).inverse().transpose() * theta2_dir;
        }
    }
}

}
