#include "RigidBodyTemplate.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

namespace bird1 {

RigidBodyTemplate::RigidBodyTemplate(Eigen::Ref<MatrixX3d> V,
                                     Eigen::Ref<MatrixX3i> F,
                                     double scale) : volume_(0), radius_(0) {
    inertiaTensor_.setZero();
    this->V = V * scale;
    this->F = F;

    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate() { }

void RigidBodyTemplate::initialize() {
    this->volume_ = computeVolume();
    this->com_ = computeCenterOfMass();

    /** Translate all vertices by the center of mass. */
    for (int v = 0; v < V.rows(); v++) V.row(v) -= this->com_;

    this->inertiaTensor_ = computeInertiaTensor();    
}

double RigidBodyTemplate::computeVolume() {
    double volume = 0;
    for (int face = 0; face < F.rows(); face++) {
        Vector3d a = V.row(F(face, 0));
        Vector3d b = V.row(F(face, 1));
        Vector3d c = V.row(F(face, 2));

        Vector3d norm_scaled = (b - a).cross(c - a);
        Vector3d norm = norm_scaled.normalized();

        /** Stoke's Theorem Magic **/
        volume += (1/6.) * norm(0) * norm_scaled.norm() * (a(0) + b(0) + c(0));
    }

    return volume;
}

Vector3d RigidBodyTemplate::computeCenterOfMass() {
    Vector3d cm(0, 0, 0);
    double volume = 0;
    for (int face = 0; face < F.rows(); face++) {
        Vector3d a = V.row(F(face, 0));
        Vector3d b = V.row(F(face, 1));
        Vector3d c = V.row(F(face, 2));

        Vector3d norm_scaled = (b - a).cross(c - a);
        Vector3d norm = norm_scaled.normalized();

        /** More Stoke's Theorem Magic. */
        volume += (1/6.) * norm(0) * norm_scaled.norm() * (a(0) + b(0) + c(0));
        cm += (1/24.) * norm_scaled.norm() * norm.cwiseProduct(
            a.cwiseProduct(a + b + c) + b.cwiseProduct(b + c) + c.cwiseProduct(c));
    }

    cm /= volume;
    return cm;
}

Matrix3d RigidBodyTemplate::computeInertiaTensor() {
    /**
     * The inertia tensor has a relatively simple form consisting of xy, xz, yz, x^2, y^2, and z^2 terms.
     * We'll use Stoke's Theorem to integrate up each term independently, and then create the final inertia tensor.
     */
    Vector3d squares = Vector3d::Zero(); // Stores x^2, y^2, z^2
    Vector3d joints = Vector3d::Zero(); // Stores yz, xz, xy

    double x3 = 0.0, y3 = 0.0, z3 = 0.0;

    for (int face = 0; face < F.rows(); face++) {
        // Need to do elementwise operations, so have array reps as well.
        Vector3d a = V.row(F(face, 0)), b = V.row(F(face, 1)), c = V.row(F(face, 2));
        Array3d ar = a.array(), br = b.array(), cr = b.array();

        Vector3d norm_scaled = (b - a).cross(c - a);
        Vector3d norm = norm_scaled.normalized();

        double a0 = a(0), a1 = a(1), a2 = a(2), b0 = b(0), b1 = b(1), b2 = b(2), c0 = c(0), c1 = c(1), c2 = c(2);

        // Using Eigen arrays to vectorize this gives slightly different results (by ~0.04) :|
        squares(0) += norm_scaled.norm() * norm(0) * (a0*a0*a0 + (a0*a0 + b0*b0 + c0*c0)*(b0+c0) + a0*(b0*b0 + b0*c0 + c0*c0)) / 60.0;
        squares(1) += norm_scaled.norm() * norm(1) * (a1*a1*a1 + (a1*a1 + b1*b1 + c1*c1)*(b1+c1) + a1*(b1*b1 + b1*c1 + c1*c1)) / 60.0;
        squares(2) += norm_scaled.norm() * norm(2) * (a2*a2*a2 + (a2*a2 + b2*b2 + c2*c2)*(b2+c2) + a2*(b2*b2 + b2*c2 + c2*c2)) / 60.0;

        // The integral of xyz over the triangle face. This one is ugly; I used the mathematica CForm command.
        double xyz = (2*a2*b0*b1 + 6*b0*b1*b2 + a2*b1*c0 + 2*b1*b2*c0 + a2*b0*c1 + 2*b0*b2*c1 + 2*a2*c0*c1 + 2*b2*c0*c1
            + 2*(b0*(b1 + c1) + c0*(b1 + 3*c1))*c2 + a0*(2*b1*b2 + b2*c1 + 2*a2*(b1 + c1) + b1*c2 + 2*c1*c2 + 2*a1*(3*a2 + b2 + c2))
            + a1*(2*a2*(b0 + c0) + b0*(2*b2 + c2) + c0*(b2 + 2*c2))) / 120.;
        joints += norm_scaled.norm() * norm * xyz;
    }

    // Now we can extract out the terms and construct the inertia tensor.
    double x2 = squares(0), y2 = squares(1), z2 = squares(2);
    double yz = joints(0), xz = joints(1), xy = joints(2);

    Matrix3d inertiaTensor;
    inertiaTensor <<
        y2 + z2, -xy, -xz,
        -xy, x2 + z2, -yz,
        -xz, -yz, x2 + y2;

    return inertiaTensor;
}

}
