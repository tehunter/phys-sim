#include "RigidBodyTemplate.h"
#include "Distance.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <map>

namespace bird2 {

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(Eigen::Ref<Eigen::MatrixX3d> V,
                                     Eigen::Ref<Eigen::MatrixX3i> F,
                                     double scale) : volume_(0)
{
    inertiaTensor_.setZero();
    Eigen::MatrixXd mV = V * scale;
    Eigen::MatrixXi mF = F;

    igl::copyleft::tetgen::tetrahedralize(mV, mF, "pq1.414a0.01", this->V, this->T, this->F);
    computeFaces();
    initialize();
}

RigidBodyTemplate::RigidBodyTemplate(const Eigen::MatrixX3d& verts, const Eigen::MatrixX4i& tets)
    : volume_(0) {
    V = verts;
    T = tets;
    computeFaces();
    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate() { }

void RigidBodyTemplate::initialize() {
    volume_ = computeVolume();
    com_ = computeCenterOfMass();
    for (int i = 0; i < V.rows(); i++)
        V.row(i) -= com_;
    inertiaTensor_ = computeInertiaTensor();
    distances_ = computeDistances();
}

void RigidBodyTemplate::computeFaces() {
    struct triple {
        triple(int aa, int bb, int cc)
            : a(aa)
            , b(bb)
            , c(cc)
        {
            if (a < b)
                std::swap(a, b);
            if (a < c)
                std::swap(a, c);
            if (b < c)
                std::swap(b, c);
        }

        int a, b, c;
        bool operator<(const triple& other) const
        {
            if (a < other.a)
                return true;
            else if (a > other.a)
                return false;
            if (b < other.b)
                return true;
            else if (b > other.b)
                return false;
            return c < other.c;
        }
    };

    int ntets = (int)T.rows();
    MatrixX3i allfaces(4 * ntets, 3);
    Matrix<int, 4, 3> faceidx;
    faceidx << 0, 1, 3,
               3, 1, 2,
               3, 2, 0,
               0, 2, 1;

    for (int i = 0; i < ntets; i++) {
        Vector4i tet = T.row(i);
        for (int face = 0; face < 4; face++) {
            for (int k = 0; k < 3; k++)
                allfaces(4 * i + face, k) = tet[faceidx(face, k)];
        }
    }

    map<triple, vector<int>> faces;
    for (int i = 0; i < 4 * ntets; i++) {
        triple t(allfaces(i, 0), allfaces(i, 1), allfaces(i, 2));
        faces[t].push_back(i / 4);
    }

    int nfaces = 0;
    for (map<triple, vector<int>>::iterator it = faces.begin(); it != faces.end(); ++it)
        if (it->second.size() == 1)
            nfaces++;

    F.resize(nfaces, 3);
    int idx = 0;

    for (int i = 0; i < 4 * ntets; i++) {
        triple t(allfaces(i, 0), allfaces(i, 1), allfaces(i, 2));
        if (faces[t].size() == 1) {
            F.row(idx) = allfaces.row(i);
            idx++;
        }
    }
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

Eigen::VectorXd RigidBodyTemplate::computeDistances() {
    int nverts = (int)V.rows();
    int nfaces = (int)F.rows();
    Eigen::VectorXd distances;
    distances.resize(nverts);
    for (int i = 0; i < nverts; i++) {
        double dist = numeric_limits<double>::infinity();
        for (int j = 0; j < nfaces; j++) {
            double dummy;
            if (Distance::vertexPlaneDistanceLessThan(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dist)) {
                Vector3d distvec = Distance::vertexFaceDistance(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dummy, dummy, dummy);
                dist = min(dist, distvec.norm());
            }
        }

        distances(i) = dist;
    }
    return distances;
}

double RigidBodyTemplate::distance(Vector3d p, int tet) const {
    // Extract the four corners of the tetrahedra that we are in.
    Vector3d v1 = V.row(T(tet, 0)), v2 = V.row(T(tet, 1)),
        v3 = V.row(T(tet, 2)), v4 = V.row(T(tet, 3));
    double d1 = distances_[T(tet, 0)], d2 = distances_[T(tet, 1)],
        d3 = distances_[T(tet, 2)], d4 = distances_[T(tet, 3)];

    // Create the T matrix so we can invert it.
    Matrix3d t = Matrix3d::Zero();
    t.col(0) = v1 - v4;
    t.col(1) = v2 - v4;
    t.col(2) = v3 - v4;

    Vector3d bary = t.inverse() * (p - v4);

    return bary(0) * d1 + bary(1) * d2 + bary(2) * d3 + (1 - bary(0) - bary(1) - bary(2)) * d4;
}

Vector3d RigidBodyTemplate::Ddistance(int tet) const {
    // Extract the four corners of the tetrahedra that we are in.
    Vector3d v1 = V.row(T(tet, 0)), v2 = V.row(T(tet, 1)),
        v3 = V.row(T(tet, 2)), v4 = V.row(T(tet, 3));
    double d1 = distances_[T(tet, 0)], d2 = distances_[T(tet, 1)],
        d3 = distances_[T(tet, 2)], d4 = distances_[T(tet, 3)];

    // Create the T matrix so we can invert it.
    Matrix3d t = Matrix3d::Zero();
    t.col(0) = v1 - v4;
    t.col(1) = v2 - v4;
    t.col(2) = v3 - v4;

    Vector3d d(d1 - d4, d2 - d4, d3 - d4);

    return d.transpose() * t.inverse();
}

}
