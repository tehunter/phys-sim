#include "GooCore.h"
#include "SimParameters.h"
#include <Eigen/Geometry>
#include <Eigen/SparseQR>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>

using namespace Eigen;

namespace goo2 {

/** A generic implementation of Newton's method which uses a sparse Jacobian. */
template<class F, class J> VectorXd newtonsSparse(VectorXd initial, int max_iters, double tolerance, const F&& func, const J&& jacobian) {
    VectorXd guess = initial;
    VectorXd value = func(initial);

    int steps = 0;
    while (steps < max_iters && value.norm() > tolerance) {
        auto jacob = jacobian(guess);
        jacob.makeCompressed();

        SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        solver.analyzePattern(jacob);
        solver.factorize(jacob);

        guess = solver.solve(-value) + guess;
        value = func(guess);
        steps++;
    }

    std::cout << "Ran for " << steps << " iterations with error " << value.norm() << std::endl;

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

/** Compute the cross product of two 2d vectors (by assuming z = 0). */
Vector3d cross2d(Vector2d a, Vector2d b) {
    Vector3d a3(a(0), a(1), 0.0);
    Vector3d b3(b(0), b(1), 0.0);

    return a3.cross(b3);
}

Vector2d cross2dWithZ(Vector2d input) {
    Vector3d input3(input(0), input(1), 0.0);
    Vector3d z(0, 0, 1.0);

    Vector3d result = input3.cross(z);
    Vector2d result2(result(0), result(1));
    return result2;
}

GooCore::GooCore() {
    params_ = std::make_shared<SimParameters>();
}

GooCore::~GooCore() { }

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> GooCore::getCurrentMesh() const {
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

    if (params_->floorEnabled) {
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
        Vector2d sourcepos = particles_[c->p1].pos;
        Vector2d destpos = particles_[c->p2].pos;
        Vector2d vec = destpos - sourcepos;
        Vector2d perp(-vec[1], vec[0]);
        perp /= perp.norm();

        double dist = (sourcepos - destpos).norm();

        Eigen::Vector3d color;
        double width;
        switch (c->getType()) {
        case SimParameters::CT_SPRING: {
            if (c->associatedBendingStencils.empty())
                color << 0.0, 0.0, 1.0;
            else
                color << 0.75, 0.5, 0.75;
            width = baselinewidth / (1.0 + 20.0 * dist * dist);

            break;
        }
        case SimParameters::CT_RIGIDROD: {
            Eigen::Vector3d color;
            if (c->associatedBendingStencils.empty())
                color << 1.0, 0.0, 1.0;
            else
                color << 1.0, 1.0, 0.3;
            width = baselinewidth;
            break;
        }
        default:
            break;
        }

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

    for (int i = 0; i < nparticles; i++) {
        double radius = baseradius * sqrt(getTotalParticleMass(i));
        radius *= (1.0 + pulsefactor * sin(pulsespeed * time_));

        Eigen::Vector3d color(0, 0, 0);

        if (particles_[i].fixed) {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++) {
            vertexColors.push_back(color);
        }

        verts.push_back(Eigen::Vector3d(particles_[i].pos[0], particles_[i].pos[1], 0));

        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++) {
            verts.emplace_back(particles_[i].pos[0] + radius * cos(2 * PI * j / numcirclewedges),
                               particles_[i].pos[1] + radius * sin(2 * PI * j / numcirclewedges),
                               0);
        }

        for (int j = 0; j <= numcirclewedges; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1)));
        }

        idx += numcirclewedges + 2;
    }

    for (const auto& saw : saws_) {
        double outerradius = saw.radius;
        double innerradius = (1.0 - sawdepth) * outerradius;

        Eigen::Vector3d color(0.5, 0.5, 0.5);

        int spokes = 2 * sawteeth;
        for (int j = 0; j < spokes + 2; j++) {
            vertexColors.push_back(color);
        }

        verts.emplace_back(saw.pos[0], saw.pos[1], 0);

        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++) {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.emplace_back(saw.pos[0] + radius * cos(2 * PI * i / spokes + sawangspeed * time_),
                               saw.pos[1] + radius * sin(2 * PI * i / spokes + sawangspeed * time_),
                               0.0);
        }

        for (int j = 0; j <= spokes; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1)));
        }

        idx += spokes + 2;
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

void GooCore::initSimulation() {
    particle_unique_id_ = 0;
    time_ = 0;
    particles_.clear();
    for (auto c: connectors_) delete c;
    connectors_.clear();
    saws_.clear();
    bendingStencils_.clear();
}

bool GooCore::simulateOneStep() {
    VectorXd config = this->configuration_vector();
    VectorXd velocity = this->velocity_vector();
    DiagonalMatrix<double, Dynamic> inv_mass = this->inv_mass_matrix();
    double h = params_->timeStep;

    // Constraint handling & time integration.
    VectorXd next_config, next_velocity;
    switch (params_->constraintHandling) {
        case SimParameters::CH_PENALTY: {
            // Penalty is automatically included in the forces, so just velocity verlet.
            next_config = config + h * velocity;
            next_velocity = velocity + h * inv_mass * this->force(next_config);
            break;
        }
        case SimParameters::CH_STEPPROJECT: {
            // Start with velocity verlet, then find nearest constrained point using Newton.
            VectorXd unconstrained_config = config + h * velocity;
            VectorXd unconstrained_velocity = velocity + h * inv_mass * this->force(unconstrained_config);

            // Correct the point using Lagrange multipliers for optimization.
            int num_constraints = this->num_rigid_rods();
            VectorXd initial_full = VectorXd::Zero(config.size() + num_constraints);
            initial_full.segment(0, config.size()) = unconstrained_config;

            VectorXd constrained_full = newtonsSparse(
                initial_full, params_->NewtonMaxIters, params_->NewtonTolerance,
                [&](const VectorXd& input) {
                    VectorXd result = VectorXd::Zero(input.size());

                    VectorXd q = input.segment(0, config.size());
                    VectorXd lambda = input.segment(config.size(), num_constraints);

                    VectorXd summation = VectorXd::Zero(config.size());
                    int index = 0;
                    for (const auto& conn : this->connectors_) {
                        if (conn->getType() != SimParameters::CT_RIGIDROD) continue;

                        RigidRod* rod = (RigidRod*) conn;
                        const Vector2d pos1 = q.segment<2>(2*rod->p1);
                        const Vector2d pos2 = q.segment<2>(2*rod->p2);

                        summation.segment<2>(2*rod->p1) += lambda(index) * 2 * (pos1 - pos2);
                        summation.segment<2>(2*rod->p2) -= lambda(index) * 2 * (pos1 - pos2);
                        result(config.size() + index) = (pos1 - pos2).squaredNorm() - rod->length * rod->length;
                        index++;
                    }

                    result.segment(0, config.size()) = (q - unconstrained_config) + inv_mass * summation;
                    return result;
                },
                [&](const VectorXd& input) {
                    VectorXd q = input.segment(0, config.size());
                    VectorXd lambda = input.segment(config.size(), num_constraints);

                    // Sum the hessian Hg_i over all rigid rods (into 'summation').
                    int index = 0;
                    SparseMatrix<double> summation(config.size(), config.size());
                    for (const auto& conn : this->connectors_) {
                        if (conn->getType() != SimParameters::CT_RIGIDROD) continue;

                        const RigidRod* rod = (RigidRod*) conn;
                        SparseMatrix<double> sel1 = this->selection_matrix(rod->p1);
                        SparseMatrix<double> sel2 = this->selection_matrix(rod->p2);
                        SparseMatrix<double> seldiff = sel1 - sel2;

                        summation += 2 * lambda(index) * seldiff.transpose() * seldiff;
                        index++;
                    }

                    SparseMatrix<double> result(input.size(), input.size());

                    // Upper left corner: compute I + M^-1 * Summation
                    // Add identity to the upper part of the result:
                    for (int i = 0; i < config.size(); i++) result.insert(i, i) = 1.0;

                    // Add inv_mass * summation to result directly.
                    SparseMatrix<double> inv_mass_summ = inv_mass * summation;
                    inv_mass_summ.resize(input.size(), input.size());
                    result += inv_mass_summ;

                    // Other two parts: copy dg to (config.size(), 0) & M^-1 * dg transpose to (0, config.size()).
                    SparseMatrix<double> dg = this->df_constraint(q);
                    for (int k = 0; k < dg.outerSize(); k++) {
                        for (SparseMatrix<double>::InnerIterator it(dg, k); it; ++it)
                            result.insert(config.size() + it.row(), it.col()) = it.value();
                    }

                    SparseMatrix<double> mdg = inv_mass * dg.transpose();
                    for (int k = 0; k < mdg.outerSize(); k++) {
                        for (SparseMatrix<double>::InnerIterator it(mdg, k); it; ++it)
                            result.insert(it.row(), config.size() + it.col()) = it.value();
                    }

                    return result;
                });
            
            next_config = constrained_full.segment(0, config.size());
            next_velocity = unconstrained_velocity + (next_config - unconstrained_config) / h;
            break;
        }
        case SimParameters::CH_LAGRANGEMULT:
            // Fancy integrator which just immediately finds a valid constrained solution.
            // Start by computing the next valid configuration; note that we are abusing 'velocity' here,
            // since this integrator actually uses momentum instead of velocity. Hopefully it still works.
            next_config = config + h * inv_mass * velocity;

            VectorXd initial = VectorXd::Zero(this->num_rigid_rods());
            VectorXd base_value = next_config + h * inv_mass * velocity + h * h * inv_mass * this->force(next_config);
            SparseMatrix<double> dg = this->df_constraint(next_config);

            VectorXd lambda = newtonsSparse(initial, params_->NewtonMaxIters, params_->NewtonTolerance,
                [&](const VectorXd& input) {
                    return this->constraint(base_value + h * h * inv_mass * dg.transpose() * input);
                },
                [&](const VectorXd& input) {
                    SparseMatrix<double> partial = h * h * inv_mass * dg.transpose();
                    SparseMatrix<double> result = this->df_constraint(base_value + partial * input) * partial;
                    return result;
                });

            // Now we compute the next lambda using the current lambda.
            next_velocity = velocity + h * this->force(next_config) + h * dg.transpose() * lambda;
            break;
    }

    // Save the timestep so it shows up in the render.
    this->persist_configuration(next_config, next_velocity);

    /** Bulk deletion of any particles, connectors, or bending stencils which go out of bounds or overlap a saw. */
    // Sets containing indices of deleted particles, connectors, and bending stencils.
    std::unordered_set<size_t> dpart, dconn, dbending;

    // Start by collecting all of the particles which are out of bounds or overlapping a saw.
    for (size_t index = 0; index < this->particles_.size(); index++) {
        Vector2d pos = this->particles_[index].pos;
        bool in_bounds = std::abs(pos[0]) <= SimParameters::SIM_DIMENSION && std::abs(pos[1]) <= SimParameters::SIM_DIMENSION;
        
        if (!in_bounds || this->overlaps_saw(pos)) dpart.insert(index);
    }

    // Now collect all connectors which have a deleted particle, or are overlapping a saw, or are overstrained.
    for (size_t index = 0; index < this->connectors_.size(); index++) {
        Connector* conn = this->connectors_[index];

        // Delete connectors connected to a deleted particle.
        if (dpart.count(conn->p1) || dpart.count(conn->p2)) {
            dconn.insert(index);
            continue;
        }

        // Handle springs which are overstrained.
        if (conn->getType() == SimParameters::ConnectorType::CT_SPRING) {
            Spring* spring = (Spring*) conn;

            double distance = (this->particles_[spring->p1].pos - this->particles_[spring->p2].pos).norm();
            double strain = (distance - spring->restlen) / spring->restlen;

            if (spring->canSnap && strain > this->params_->maxSpringStrain) {
                dconn.insert(index);
                continue;
            }
        }

        // Finally, check if connector overlaps a saw, if so then delete it.
        if (this->overlaps_saw(this->particles_[conn->p1].pos, this->particles_[conn->p2].pos)) dconn.insert(index);
    }

    // And finally collect all bending stencils connected to a deleted particle or spring.
    for (size_t index : dconn) {
        for (int bending : this->connectors_[index]->associatedBendingStencils) dbending.insert(bending);
    }

    // Early return if we don't need to delete anything.
    if (dpart.size() == 0 && dconn.size() == 0 && dbending.size() == 0) return false;
    // TODO: Simpler logic if only connectors have been deleted, since reindexing is not required.

    // Vectors which map old particle indices -> new particle indices, and same for bending stencils.
    std::vector<size_t> ipart(this->particles_.size()), ibending(this->bendingStencils_.size());
    size_t valid_index = 0;
    for (size_t index = 0; index < this->particles_.size(); index++) {
        if (dpart.count(index)) ipart[index] = -1;
        else ipart[index] = valid_index++;
    }

    valid_index = 0;
    for (size_t index = 0; index < this->bendingStencils_.size(); index++) {
        if (dbending.count(index)) ibending[index] = -1;
        else ibending[index] = valid_index++;
    }

    // Now construct new vectors of particles, connectors, and bending stencils.
    std::vector<Particle, Eigen::aligned_allocator<Particle>> new_particles;
    new_particles.reserve(this->particles_.size() - dpart.size());
    for (size_t index = 0; index < this->particles_.size(); index++) {
        if (ipart[index] == -1) continue;
        new_particles.push_back(this->particles_[index]);
    }

    std::vector<Connector*> new_connectors;
    new_connectors.reserve(this->connectors_.size() - dconn.size());
    for (size_t index = 0; index < this->connectors_.size(); index++) {
        if (dconn.count(index)) continue;

        Connector* conn = this->connectors_[index];
        // Update connector particle indexes, as well as bending stencil indexes.
        conn->p1 = ipart[conn->p1];
        conn->p2 = ipart[conn->p2];

        std::set<int> new_stencils;
        for (int old_bending : conn->associatedBendingStencils) {
            if (ibending[old_bending] == -1) continue;
            new_stencils.insert(ibending[old_bending]);
        }
        conn->associatedBendingStencils = std::move(new_stencils);

        new_connectors.push_back(conn);
    }

    std::vector<BendingStencil> new_bending;
    new_bending.reserve(this->bendingStencils_.size() - dbending.size());
    for (size_t index = 0; index < this->bendingStencils_.size(); index++) {
        if (ibending[index] == -1) continue;
        BendingStencil stencil = this->bendingStencils_[index];
        stencil.p1 = ipart[stencil.p1];
        stencil.p2 = ipart[stencil.p2];
        stencil.p3 = ipart[stencil.p3];

        new_bending.push_back(stencil);
    }

    this->particles_ = std::move(new_particles);
    this->connectors_ = std::move(new_connectors);
    this->bendingStencils_ = std::move(new_bending);
    return false;
}

VectorXi GooCore::addParticle(double x, double y) {
    int new_index = particles_.size();
    int new_id = particle_unique_id_++;
    Vector2d new_position(x,y);

    double mass = params_->particleFixed ? std::numeric_limits<double>::infinity() : params_->particleMass;
    particles_.emplace_back(new_position, mass, params_->particleFixed, false, new_id);

    // Elastic rods create multiple particles, so we collect all of their indexes and return them here.
    std::vector<int> ret_ids;
    ret_ids.push_back(new_id);

    // Iterate through all other particles to find particles within the max spring dist.
    for(int index = 0; index < particles_.size(); index++) {
        if (index == new_index) continue;
        if (particles_[index].inert) continue;

        Vector2d target_position = particles_[index].pos;
        double distance = (target_position - new_position).norm();
        if (distance >= params_->maxSpringDist) continue;

        // TODO: Spring mass and can snap?
        switch (params_->connectorType) {
            case SimParameters::ConnectorType::CT_SPRING:
                connectors_.push_back(new Spring(new_index, index, 0.0, params_->springStiffness / distance, distance, true));
                break;
            case SimParameters::ConnectorType::CT_RIGIDROD:
                connectors_.push_back(new RigidRod(new_index, index, 0.0, distance));
                break;
            default: {
                // Create <num_segments> - 1 intermediate particles.
                std::vector<Vector2d> positions;
                positions.reserve(params_->rodSegments + 1);
                std::vector<int> indexes;
                indexes.reserve(params_->rodSegments + 1);
                positions.push_back(new_position);
                indexes.push_back(new_index);

                for (int np = 0; np < params_->rodSegments - 1; np++) {
                    int partial_new_id = particle_unique_id_++;
                    int partial_new_index = particles_.size();
                    Vector2d partial_position = new_position + (target_position - new_position) / params_->rodSegments * (np+1);
                    particles_.emplace_back(partial_position, 0.0, false, true, partial_new_id);

                    // Add this particle to the return, as well as intermediate tracking stuff.
                    ret_ids.push_back(partial_new_id);
                    indexes.push_back(partial_new_index);
                    positions.push_back(partial_position);
                }

                indexes.push_back(index);
                positions.push_back(target_position);

                // Add connectors to pairwise elements, and bending stencils to the associated connectors.
                std::vector<Connector*> created_conns;
                for (int i = 0; i < indexes.size() - 1; i++) {
                    double distance = (positions[i] - positions[i + 1]).norm();
                    connectors_.push_back(new Spring(indexes[i], indexes[i+1], distance * params_->rodDensity,
                        params_->rodStretchingStiffness / distance, distance, false));
                    created_conns.push_back(connectors_[connectors_.size() - 1]);
                }

                // Add bending stencils to triples of elements.
                for (int i = 0; i < indexes.size() - 2; i++) {
                    double d1 = (positions[i] - positions[i + 1]).norm();
                    double d2 = (positions[i + 1] - positions[i + 2]).norm();

                    bendingStencils_.emplace_back(indexes[i], indexes[i+1], indexes[i+2], 2*params_->rodBendingStiffness / (d1 + d2));
                    created_conns[i]->associatedBendingStencils.insert(bendingStencils_.size()-1);
                    created_conns[i+1]->associatedBendingStencils.insert(bendingStencils_.size()-1);
                }

                break;
            }
        }
    }

    VectorXi result(ret_ids.size());
    for (int i = 0; i < ret_ids.size(); i++) result(i) = ret_ids[i];
    return result;
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

VectorXi GooCore::queryConnectivity(VectorXi from, VectorXi to) {
    VectorXi result = VectorXi::Zero(from.size());
    std::unordered_map<int, int> uid2index;
    for (size_t i = 0; i < particles_.size(); i++)
        uid2index[particles_[i].uid] = static_cast<int>(i);
    
    for (int i = 0; i < from.size(); i++) {
        int findex = uid2index[from(i)];
        int tindex = uid2index[to(i)];

        for (const auto& conn : connectors_) {
            if ((conn->p1 == findex && conn->p2 == tindex)
                || (conn->p2 == findex && conn->p1 == tindex)) {
                result(i) = 1;
                break;
            }
        }
    }

    return result;
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
    VectorXd masses = this->particle_masses();
    for (size_t i = 0; i < particles_.size(); i++) {
        matrix.diagonal()[2*i] = masses(i);
        matrix.diagonal()[2*i+1] = masses(i);
    }

    return matrix;
}

DiagonalMatrix<double, Dynamic> GooCore::inv_mass_matrix() const {
    // The inverse of a diagonal matrix is the elementwise reciprocal.
    VectorXd masses = this->particle_masses();
    DiagonalMatrix<double, Dynamic> matrix(2*particles_.size());
    for (size_t i = 0; i < particles_.size(); i++) {
        if (particles_[i].fixed) {
            matrix.diagonal()[2*i] = 0.0;
            matrix.diagonal()[2*i+1] = 0.0;
        } else {
            matrix.diagonal()[2*i] = 1.0/masses(i);
            matrix.diagonal()[2*i+1] = 1.0/masses(i);
        }
    }

    return matrix;
}

VectorXd GooCore::force_gravity(const VectorXd& config) const {
    // TODO: Consider passing in the particle masses seperately.
    // Gravity is a constant force downward in the Y direction.
    VectorXd masses = this->particle_masses();
    VectorXd result = VectorXd::Zero(config.size());
    for (int i = 1; i < config.size(); i += 2) {
        if (particles_[(i - 1)/2].fixed) continue;
        result(i) = params_->gravityG * masses((i - 1)/2);
    }

    return result;
}

/** Compute the force vector for the springs. */
VectorXd GooCore::force_spring(const VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (auto conn : this->connectors_) {
        if (conn->getType() != SimParameters::ConnectorType::CT_SPRING) continue;

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

VectorXd GooCore::force_floor(const VectorXd& config) const {
    // A simple quadratic potential/linear force which dissuades particles from going below the floor too long.
    VectorXd result = VectorXd::Zero(config.size());
    VectorXd masses = this->particle_masses();
    for (int i = 1; i < config.size(); i += 2) {
        if (particles_[(i - 1)/2].fixed) continue;
        double pradius = 0.02 * sqrt(masses((i - 1)/2));
        double pbottom = config[i] - pradius;
        if (pbottom > -0.5) continue;

        double offset = (-0.50 - pbottom);
        result[i] = SimParameters::FLOOR_STRENGTH * offset;
    }
    
    return result;
}

VectorXd GooCore::force_damping(const Eigen::VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (auto conn : this->connectors_) {
        if (conn->getType() != SimParameters::ConnectorType::CT_SPRING) continue;

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

/** Compute the penalty force (and derivative) for rigid rods. */
VectorXd GooCore::force_rigid_penalty(const VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (auto conn : this->connectors_) {
        if (conn->getType() != SimParameters::ConnectorType::CT_RIGIDROD) continue;

        const RigidRod* rod = (RigidRod*) conn;
        const Vector2d pos1 = config.segment<2>(2*rod->p1);
        const Vector2d pos2 = config.segment<2>(2*rod->p2);

        double constraint = (pos1 - pos2).squaredNorm() - rod->length * rod->length;
        Vector2d penalty = -4 * params_->penaltyStiffness * constraint * (pos1 - pos2);

        result.segment<2>(2*rod->p1) += penalty;
        result.segment<2>(2*rod->p2) -= penalty;
    }

    return result;
}

/** Compute the elastic bending force */
VectorXd GooCore::force_bending(const VectorXd& config) const {
    VectorXd result = VectorXd::Zero(config.size());
    for (const auto& bending : this->bendingStencils_) {
        const Vector2d pi = config.segment<2>(2*bending.p1);
        const Vector2d pj = config.segment<2>(2*bending.p2);
        const Vector2d pk = config.segment<2>(2*bending.p3);

        Vector2d u = (pj - pi), v = (pk - pj);

        // u = p_j - p_i, v = p_k - p_j
        // theta = 2atan2((u cross v) * z_hat, norm(u)norm(v)+(u dot v))
        float theta = 2 * std::atan2(cross2d(u,v)[2], u.norm()*v.norm() + u.dot(v));

        Vector2d fi = bending.kb * theta * cross2dWithZ(u).transpose() / u.squaredNorm();
        Vector2d fk = bending.kb * theta * cross2dWithZ(v).transpose() / v.squaredNorm();
        Vector2d fj = -fi - fk;

        result.segment<2>(2*bending.p1) += fi;
        result.segment<2>(2*bending.p2) += fj;
        result.segment<2>(2*bending.p3) += fk;
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
    if (params_->constraintHandling == SimParameters::CH_PENALTY)
        result += this->force_rigid_penalty(config);
    if (params_->bendingEnabled) result += this->force_bending(config);

    // Clear forces on fixed particles.
    for (int i = 0; i < particles_.size(); i++) {
        if (particles_[i].fixed) result.segment<2>(2*i) = Vector2d::Zero();
    }

    return result;
}

/** Compute the value of the constraint function g(q). */
VectorXd GooCore::constraint(const Eigen::VectorXd& config) const {
    VectorXd result = VectorXd::Zero(this->num_rigid_rods());
    int index = 0;
    for (const auto& conn : this->connectors_) {
        if (conn->getType() != SimParameters::CT_RIGIDROD) continue;

        RigidRod* rod = (RigidRod*) conn;
        const Vector2d pos1 = config.segment<2>(2*rod->p1);
        const Vector2d pos2 = config.segment<2>(2*rod->p2);

        result(index) = (pos1 - pos2).squaredNorm() - rod->length * rod->length;
        index++;
    }

    return result;
}

/** Compute the derivative of the constraint function g(q). */
SparseMatrix<double> GooCore::df_constraint(const Eigen::VectorXd& config) const {
    int num_constraints = this->num_rigid_rods();
    SparseMatrix<double> dg(num_constraints, config.size());
    int index = 0;
    for (const auto& conn : this->connectors_) {
        if (conn->getType() != SimParameters::CT_RIGIDROD) continue;

        RigidRod* rod = (RigidRod*) conn;
        const Vector2d pos1 = config.segment<2>(2*rod->p1);
        const Vector2d pos2 = config.segment<2>(2*rod->p2);
        const Vector2d quant = 2 * (pos1 - pos2);

        /**
        * The sparse matrix equivalent of
        *  dg.row(index).segment<2>(2*rod->p1) += 2 * (pos1 - pos2);
        *  dg.row(index).segment<2>(2*rod->p2) -= 2 * (pos1 - pos2);
        */

        dg.coeffRef(index, 2*rod->p1) += quant[0];
        dg.coeffRef(index, 2*rod->p1+1) += quant[1];
        dg.coeffRef(index, 2*rod->p2) -= quant[0];
        dg.coeffRef(index, 2*rod->p2+1) -= quant[1];

        index++;
    }

    return dg;
}

/** Utility method for creating a selection matrix so I don't go insane. */
SparseMatrix<double> GooCore::selection_matrix(int index) const {
    SparseMatrix<double> result(2, 2*this->particles_.size());
    result.insert(0, 2*index) = 1.0;
    result.insert(1, 2*index+1) = 1.0;
    return result;
}

/** Compute the masses for all particles efficiently in one iteration over the connectors. */
VectorXd GooCore::particle_masses() const {
    VectorXd result(particles_.size());
    for (int index = 0; index < particles_.size(); index++)
        result(index) = particles_[index].mass;
    
    for (auto& conn : connectors_) {
        result(conn->p1) += conn->mass/2;
        result(conn->p2) += conn->mass/2;
    }

    return result;
}

/** Return true if the line segment overlaps a saw, and false otherwise. */
bool GooCore::overlaps_saw(Vector2d pos1, Vector2d pos2) const {
    for (size_t si = 0; si < this->saws_.size(); si++) {
        if (segment_distance(pos1, pos2, this->saws_[si].pos) <= this->saws_[si].radius) return true;
    }

    return false;
}

/** Return true if the given point overlaps a saw, and false otherwise. */
bool GooCore::overlaps_saw(Vector2d point) const {
    for (size_t si = 0; si < this->saws_.size(); si++) {
        if ((this->saws_[si].pos - point).norm() < this->saws_[si].radius) return true;
    }

    return false;
}

}
