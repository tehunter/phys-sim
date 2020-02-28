#include "GooCore.h" 
#include <stdint.h>
#include <iostream>

#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

namespace goo1 {

/** A generic implementation of Newton's method. */
template<class F, class J> VectorXd newtons(VectorXd initial, int max_iters, double tolerance, const F&& func, const J&& jacobian) {
    VectorXd guess = initial;
    VectorXd value = func(initial);

    int steps = 0;
    while (steps < max_iters && value.norm() > tolerance) {
        auto jacob = jacobian(guess);

        guess = jacob.llt().solve(-value) + guess;
        value = func(guess);
        steps++;
    }

    return guess;
}

/** Helper function which computes the distance from the line segment (vw) to the point p. */
double segment_distance(Vector2d v, Vector2d w, Vector2d p) {
    // Parameterize the line segment by c(t) = w + t(v - w)
    const double segment_length = (v - w).norm();
    double t = (p - w).dot(v - w) / (segment_length * segment_length);
    t = std::max(0.0, std::min(t, 1.0));

    const Vector2d point = w + t * (v - w);
    return (point - p).norm();
}

GooCore::GooCore() {
	params_ = std::make_shared<SimParameters>();
}

GooCore::~GooCore() { }

void GooCore::initSimulation() {
    particle_unique_id_ = 0;
    time_ = 0;
    particles_.clear();
    for (auto c: connectors_) delete c;
    connectors_.clear();
    saws_.clear();
}

bool GooCore::simulateOneStep() {
    VectorXd config = this->configuration_vector();
    VectorXd velocity = this->velocity_vector();
    DiagonalMatrix<double, Dynamic> mass = this->mass_matrix();
    DiagonalMatrix<double, Dynamic> inv_mass = this->inv_mass_matrix();
    double h = params_->timeStep;

    // Switch on the integrator to obtain the next timestep.
    VectorXd next_config, next_velocity;
    switch (params_->integrator) {
        case SimParameters::TimeIntegrator::TI_EXPLICIT_EULER:
            next_config = config + h * velocity;
            next_velocity = velocity + h * inv_mass * this->force(config);
            break;
        case SimParameters::TimeIntegrator::TI_IMPLICIT_EULER:
            next_config = newtons(config, params_->NewtonMaxIters, params_->NewtonTolerance,
                [&](const VectorXd& vec) {
                    VectorXd result = mass * (vec - config - h * velocity);
                    result -= h * h * this->force(vec);
                    return result;
                },
                [&](const VectorXd& vec) {
                    MatrixXd result = MatrixXd(mass);
                    result -= h * h * this->dforce(vec).transpose();
                    return result;
                });

            next_velocity = velocity + h * inv_mass * this->force(next_config);
            break;
        case SimParameters::TimeIntegrator::TI_IMPLICIT_MIDPOINT:
            next_config = newtons(config, params_->NewtonMaxIters, params_->NewtonTolerance,
                [&](const VectorXd& vec) {
                    VectorXd result = mass * (vec - config - h * velocity);
                    result -= (h * h) / 2 * this->force((vec + config)/2);
                    return result;
                },
                [&](const VectorXd& vec) {
                    MatrixXd result = MatrixXd(mass);
                    result -= (h * h)/4 * this->dforce((vec + config)/2).transpose();
                    return result;
                });

            next_velocity = velocity + h * inv_mass * this->force((next_config + config)/2);
            break;
        case SimParameters::TimeIntegrator::TI_VELOCITY_VERLET:
            next_config = config + h * velocity;
            next_velocity = velocity + h * inv_mass * this->force(next_config);
            break;
    }

    // Save the timestep so it shows up in the render.
    this->persist_configuration(next_config, next_velocity);

    // Now delete any springs which are too strained or are overlapping a saw.
    for (ssize_t index = this->connectors_.size() - 1; index >= 0; index--) {
        Spring* spring = (Spring*) this->connectors_[index];
        double distance = (this->particles_[spring->p1].pos - this->particles_[spring->p2].pos).norm();
        double strain = (distance - spring->restlen) / spring->restlen;

        // If the spring is strained, delete it.
        if (strain > this->params_->maxSpringStrain) {
            this->connectors_.erase(this->connectors_.begin() + index);
            continue;
        }

        const Vector2d pos1 = this->particles_[spring->p1].pos;
        const Vector2d pos2 = this->particles_[spring->p2].pos;

        // Check if spring overlaps a saw, if so then delete it.
        for (ssize_t si = this->saws_.size() - 1; si >= 0; si--) {
            if (segment_distance(pos1, pos2, this->saws_[si].pos) <= this->saws_[si].radius) {
                this->connectors_.erase(this->connectors_.begin() + index);
                break;
            }
        }
    }

    // Delete particles which are overlapping a saw or which go out of bounds.
    for (ssize_t index = this->particles_.size() - 1; index >= 0; index--) {
        Vector2d pos = this->particles_[index].pos;
        double radius = 0.02 * sqrt(this->getTotalParticleMass(index));
        bool in_bounds = std::abs(pos[0]) <= SimParameters::SIM_DIMENSION && std::abs(pos[1]) <= SimParameters::SIM_DIMENSION;
        bool in_saw = false;
        for (ssize_t si = this->saws_.size() - 1; si >= 0; si--) {
            if ((this->saws_[si].pos - pos).norm() < this->saws_[si].radius) {
                in_saw = true;
                break;
            }
        }

        // This particle is safe, for now...
        if (in_bounds && !in_saw) continue;

        // Delete all springs connected to this particle; shift the indexes of other springs.
        this->particles_.erase(this->particles_.begin() + index);
        for (ssize_t ci = this->connectors_.size() - 1; ci >= 0; ci--) {
            Spring* spring = (Spring*) this->connectors_[ci];

            // Delete springs connected to the particle.
            if (spring->p1 == index || spring->p2 == index) {
                this->connectors_.erase(this->connectors_.begin() + ci);
                continue;
            }

            // Shift the indexes of springs not connected to the particle by one.
            if (spring->p1 >= index) spring->p1--;
            if (spring->p2 >= index) spring->p2--;
        }
    }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> GooCore::getCurrentMesh() const {
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

    double baselinewidth = 0.005;

    int numcirclewedges = 20;

    // this is terrible. But, easiest to get up and running

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    int idx = 0;

    double eps = 1e-4;


    if(params_->floorEnabled) {
        for (int i = 0; i < 6; i++) {
            vertexColors.emplace_back(0.3, 1.0, 0.3);
        }

        verts.emplace_back(-1, -0.5, eps);
        verts.emplace_back(1, -0.5, eps);
        verts.emplace_back(-1, -1, eps);

        faces.emplace_back(idx, idx + 1, idx + 2);

        verts.emplace_back(-1, -1, eps);
        verts.emplace_back(1, -0.5, eps);
        verts.emplace_back(1, -1, eps);
        faces.emplace_back(idx + 3, idx + 4, idx + 5);
        idx += 6;
    }

    for (auto c : connectors_) {
        Eigen::Vector3d color;
        if (c->associatedBendingStencils.empty())
            color << 0.0, 0.0, 1.0;
        else
            color << 0.75, 0.5, 0.75;
        Vector2d sourcepos = particles_[c->p1].pos;
        Vector2d destpos = particles_[c->p2].pos;

        Vector2d vec = destpos - sourcepos;
        Vector2d perp(-vec[1], vec[0]);
        perp /= perp.norm();

        double dist = (sourcepos - destpos).norm();

        double width = baselinewidth / (1.0 + 20.0 * dist * dist);

        for (int i = 0; i < 4; i++)
            vertexColors.emplace_back(color);

        verts.emplace_back(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps);
        verts.emplace_back(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps);
        verts.emplace_back(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps);
        verts.emplace_back(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps);

        faces.emplace_back(idx, idx + 1, idx + 2);
        faces.emplace_back(idx + 2, idx + 1, idx + 3);
        idx += 4;
    }

    int nparticles = particles_.size();

    for(int i=0; i<nparticles; i++) {
        double radius = baseradius*sqrt(getTotalParticleMass(i));
        radius *= (1.0 + pulsefactor*sin(pulsespeed*time_));

        Eigen::Vector3d color(0,0,0);

        if(particles_[i].fixed) {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++) {
            vertexColors.push_back(color);
        }


        verts.push_back(Eigen::Vector3d(particles_[i].pos[0], particles_[i].pos[1], 0));

        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++) {
            verts.emplace_back(particles_[i].pos[0] + radius * cos(2 * PI*j / numcirclewedges),
                               particles_[i].pos[1] + radius * sin(2 * PI*j / numcirclewedges),
			       0);
        }

        for (int j = 0; j <= numcirclewedges; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1)));
        }

        idx += numcirclewedges + 2;
    }

    for(const auto& saw : saws_) {
        double outerradius = saw.radius;
        double innerradius = (1.0-sawdepth)*outerradius;

        Eigen::Vector3d color(0.5,0.5,0.5);

        int spokes = 2*sawteeth;
        for (int j = 0; j < spokes + 2; j++)
        {
            vertexColors.push_back(color);
        }

        verts.emplace_back(saw.pos[0], saw.pos[1], 0);

        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++) {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.emplace_back(saw.pos[0] + radius * cos(2 * PI*i / spokes + sawangspeed*time_),
                               saw.pos[1] + radius * sin(2 * PI*i / spokes + sawangspeed*time_),
                               0.0);
        }

        for (int j = 0; j <= spokes; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1)));
        }

        idx += spokes + 2;
    }

    renderQ.resize(verts.size(),3);
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

uint32_t GooCore::addParticle(double x, double y) {
    Vector2d newpos(x,y);
    double mass = params_->particleMass;

    int newid = particles_.size();
    int ret = particle_unique_id_++;
    particles_.emplace_back(newpos, mass, params_->particleFixed, false, ret);

    // Iterate through all other particles to find particles within the max spring dist.
    for(int index = 0; index < particles_.size(); index++) {
        if (index == newid) continue;
        double distance = (particles_[index].pos - newpos).norm();
        if (distance >= params_->maxSpringDist) continue;

        // TODO: Spring mass and can snap?
        connectors_.push_back(new Spring(newid, index, 0.0, params_->springStiffness / distance, distance, true));
    }

    return ret;
}

void GooCore::addSaw(double x, double y) {
    Vector2d position(x, y);
    this->saws_.emplace_back(position, this->params_->sawRadius);
}

double GooCore::getTotalParticleMass(int idx) const {
    double mass = particles_[idx].mass;
    for(const auto c : connectors_) {
        if(c->p1 == idx || c->p2 == idx)
            mass += 0.5 * c->mass;
    }
    return mass;
}

// Obtain the configuration vector for the current state. Ordered the same as particles_ is.
VectorXd GooCore::configuration_vector() const {
    VectorXd result(2 * particles_.size());
    for (size_t index = 0; index < particles_.size(); index++) {
        result[2*index] = particles_[index].pos[0];
        result[2*index+1] = particles_[index].pos[1];
    }

    return result;
}

// Obtain the velocity vector for the current state. Ordered the same as particles_ is.
VectorXd GooCore::velocity_vector() const {
    VectorXd result(2 * particles_.size());
    for (size_t index = 0; index < particles_.size(); index++) {
        result[2*index] = particles_[index].vel[0];
        result[2*index+1] = particles_[index].vel[1];
    }

    return result;
}

// Copy position and velocity information back to the particles data structure
// so it can be rendered.
void GooCore::persist_configuration(const VectorXd& config, const VectorXd& velocity) {
    for (size_t index = 0; index < particles_.size(); index++) {
        particles_[index].prevpos = particles_[index].pos;
        particles_[index].pos = Vector2d(config[2*index], config[2*index+1]);
        particles_[index].vel = Vector2d(velocity[2*index], velocity[2*index+1]);
    }
}

DiagonalMatrix<double, Dynamic> GooCore::mass_matrix() const {
    // Just place masses along the diagonal.
    DiagonalMatrix<double, Dynamic> matrix(2*particles_.size());
    for (size_t i = 0; i < particles_.size(); i++) {
        matrix.diagonal()[2*i] = particles_[i].mass;
        matrix.diagonal()[2*i+1] = particles_[i].mass;
    }

    return matrix;
}

DiagonalMatrix<double, Dynamic> GooCore::inv_mass_matrix() const {
    // The inverse of a diagonal matrix is the elementwise reciprocal.
    DiagonalMatrix<double, Dynamic> matrix(2*particles_.size());
    for (size_t i = 0; i < particles_.size(); i++) {
        matrix.diagonal()[2*i] = 1.0/particles_[i].mass;
        matrix.diagonal()[2*i+1] = 1.0/particles_[i].mass;
    }

    return matrix;
}

VectorXd GooCore::force_gravity(const VectorXd& config) const {
    // TODO: Consider passing in the particle masses seperately.
    // Gravity is a constant force downward in the Y direction.
    VectorXd result = VectorXd::Zero(config.size());
    for (int i = 1; i < config.size(); i += 2) {
        result[i] = params_->gravityG * particles_[(i - 1)/2].mass;
    }
    
    return result;
}

MatrixXd GooCore::df_gravity(const VectorXd& config) const {
    // Amusingly, the second derivative of gravity potential energy is 0...
    return MatrixXd::Zero(config.size(), config.size());
}

/** Compute the force vector for the springs. */
VectorXd GooCore::force_spring(const VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (auto conn : this->connectors_) {
        // Move this computation to Spring class later on, and stop casting to Spring always.
        const Spring* spring = (Spring*) conn;
        const Vector2d pos1 = config.segment<2>(2*spring->p1);
        const Vector2d pos2 = config.segment<2>(2*spring->p2);

        double dist = std::max(1e-6, (pos1 - pos2).norm());
        double coeff = -spring->stiffness * (dist - spring->restlen) / dist;
        Vector2d force = coeff * (pos1 - pos2);

        result.segment<2>(2*spring->p1) += force;
        result.segment<2>(2*spring->p2) -= force;
    }

    return result;
}

/** compute the derivative of the force vector for springs. */
MatrixXd GooCore::df_spring(const VectorXd& config) const {
    // This is going to be an ugly one...
    MatrixXd result = MatrixXd::Zero(config.size(), config.size());
    for (auto conn : this->connectors_) {
        const Spring* spring = (Spring*) conn;
        const Vector2d pos1 = config.segment<2>(2*spring->p1);
        const Vector2d pos2 = config.segment<2>(2*spring->p2);
        const Vector2d diff = pos1 - pos2;
        double dist = std::max(1e-6, diff.norm());
        const Vector2d hdiff = diff / dist;

        // Explicitly using selection matrices since the math makes me sad otherwise. :(
        MatrixXd p1s = this->selection_matrix(spring->p1);
        MatrixXd p2s = this->selection_matrix(spring->p2);

        result -= spring->stiffness * (dist - spring->restlen) / dist * (p1s - p2s).transpose() * (p1s - p2s);
        result -= spring->stiffness * (p1s - p2s).transpose() * (pos1 - pos2)
            * (hdiff.transpose() / dist - hdiff.transpose() * (dist - spring->restlen) / (dist * dist)) * (p1s - p2s);
    }

    return result;
}

VectorXd GooCore::force_floor(const VectorXd& config) const {
    // A simple quadratic potential/linear force which dissuades particles from going below the floor too long.
    VectorXd result = VectorXd::Zero(config.size());
    for (int i = 1; i < config.size(); i += 2) {
        if (particles_[(i - 1)/2].fixed) continue;
        double pradius = 0.02 * sqrt(this->getTotalParticleMass((i-1)/2));
        double pbottom = config[i] - pradius;
        if (pbottom > -0.5) continue;

        double offset = (-0.50 - pbottom);
        result[i] = SimParameters::FLOOR_STRENGTH * offset;
    }
    
    return result;
}

MatrixXd GooCore::df_floor(const VectorXd& config) const {
    MatrixXd result = MatrixXd::Zero(config.size(), config.size());
    for (int i = 1; i < config.size(); i += 2) {
        if (particles_[(i - 1)/2].fixed) continue;
        double pradius = 0.02 * sqrt(this->getTotalParticleMass((i-1)/2));
        double pbottom = config[i] - pradius;
        if (pbottom > -0.5) continue;

        result(i, i) = -SimParameters::FLOOR_STRENGTH;
    }
    
    return result;
}

VectorXd GooCore::force_damping(const Eigen::VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (auto conn : this->connectors_) {
        const Spring* spring = (Spring*) conn;
        const Vector2d pos1 = config.segment<2>(2*spring->p1);
        const Vector2d pos2 = config.segment<2>(2*spring->p2);
        const Vector2d prevpos1 = this->particles_[spring->p1].prevpos;
        const Vector2d prevpos2 = this->particles_[spring->p2].prevpos;

        double k = this->params_->dampingStiffness, h = this->params_->timeStep;
        Vector2d damping = k/h * ((pos2 - prevpos2) - (pos1 - prevpos1)).transpose();

        result.segment<2>(2*spring->p1) += damping;
        result.segment<2>(2*spring->p2) -= damping;
    }

    return result;
}

MatrixXd GooCore::df_damping(const Eigen::VectorXd& config) const {
    MatrixXd result = MatrixXd::Zero(config.size(), config.size());
    for (auto conn : this->connectors_) {
        const Spring* spring = (Spring*) conn;

        double k = this->params_->dampingStiffness, h = this->params_->timeStep;

        MatrixXd sel1 = this->selection_matrix(spring->p1);
        MatrixXd sel2 = this->selection_matrix(spring->p2);

        result += -k/h * (sel1 - sel2).transpose() * (sel1 - sel2);
    }

    return result;
}

/** Computes the force in the given configuration (according to the enabled forces). */
VectorXd GooCore::force(const VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    if (params_->gravityEnabled) result += this->force_gravity(config);
    if (params_->springsEnabled) result += this->force_spring(config);
    if (params_->floorEnabled) result += this->force_floor(config);
    if (params_->dampingEnabled) result += this->force_damping(config);

    // Clear forces on fixed particles.
    for (int i = 0; i < particles_.size(); i++) {
        if (particles_[i].fixed) result.segment<2>(2*i) = Vector2d::Zero();
    }

    return result;
}

/** Computes the derivative of force in the given configuration (according to enabled forces). */
MatrixXd GooCore::dforce(const VectorXd& config) const {
    MatrixXd result = MatrixXd::Zero(config.size(), config.size());
    if (params_->gravityEnabled) result += this->df_gravity(config);
    if (params_->springsEnabled) result += this->df_spring(config);
    if (params_->floorEnabled) result += this->df_floor(config);
    if (params_->dampingEnabled) result += this->df_damping(config);

    // Clear forces on fixed particles.
    for (int i = 0; i < particles_.size(); i++) {
        if (particles_[i].fixed) {
            result.block(2*i, 0, 2, config.size()).setZero();
            result.block(0, 2*i, config.size(), 2).setZero();
        }
    }

    return result;
}

/** Utility method for creating a selection matrix so I don't go insane. */
MatrixXd GooCore::selection_matrix(int index) const {
    MatrixXd result = MatrixXd::Zero(2, 2*this->particles_.size());
    result.block<2, 2>(0, 2*index) = Matrix2d::Identity();
    return result;
}

}
