#ifndef PSYM_CORE_GOO2_GOOCORE_H
#define PSYM_CORE_GOO2_GOOCORE_H

#include "../PhysicsCore.h"
#include "SceneObjects.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

#include <memory>
#include <stdint.h>

namespace goo2 {

class SimParameters;

class GooCore : public PhysicsCore {
public:
    GooCore();
    ~GooCore();

    virtual void initSimulation() override;
    virtual bool simulateOneStep() override;

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> getCurrentMesh() const override;
    std::shared_ptr<SimParameters> getPointerToSimParameters() { return params_; }
    SimParameters getSimParameters() const { return *params_; }
    void setSimParameters(const SimParameters& nsimp) { *params_ = nsimp; }

    /*
     * queryParticle: locate the particle with specific UID
     *
     * Returns: Tuple (particle, valid).
     *          valid must be false if the particle with given uid cannot be
     *          found.
     */
    std::tuple<Particle, bool> queryParticle(uint32_t uid) const {
        for(const auto& p : particles_)
            if (p.uid == uid) return std::make_tuple(p, true);

        return std::make_tuple(Particle(Eigen::Vector2d(0.0,0.0), -1, false, false, -1), false);
    }

    /**
     * queryConnectivity: Detect whether from(i) is connected with to(i).
     * 
     * Returns: Eigen::VectorXi; ret(i) is nonzero is from(i) is connected to to(i).
     */
    Eigen::VectorXi queryConnectivity(Eigen::VectorXi from, Eigen::VectorXi to);

    /*
     * addParticle: add a particle at (x, y)
     *
     * Returns: a list of non-negative IDs which uniquely identify  the newly
     *          added particles during the lifetime of the simulation (i.e.
     *          until initSimulation() is called)
     *
     * Note: A flex rod returns mutliple particles.
     */
    Eigen::VectorXi addParticle(double x, double y);

    /* Add a saw at (x, y) */
    void addSaw(double x, double y);

    /** Obtain the configuration vector (containing particle X/Y positions). */
    Eigen::VectorXd configuration_vector() const;

    /** Obtain the configurational velocity vector (parallel to the config vector). */
    Eigen::VectorXd velocity_vector() const;

    /** Persist the data in a configuration vector. */
    void persist_configuration(const Eigen::VectorXd& config, const Eigen::VectorXd& velocity);

    /** Compute the force vector (and derivative) for gravity. */
    Eigen::VectorXd force_gravity(const Eigen::VectorXd& config) const;
    Eigen::MatrixXd df_gravity(const Eigen::VectorXd& config) const;

    /** Compute the force vector (and derivative) for the springs. */
    Eigen::VectorXd force_spring(const Eigen::VectorXd& config) const;
    Eigen::MatrixXd df_spring(const Eigen::VectorXd& config) const;

    /** Compute the force vector and derivative vector for the floor. */
    Eigen::VectorXd force_floor(const Eigen::VectorXd& config) const;
    Eigen::MatrixXd df_floor(const Eigen::VectorXd& config) const;

    /** Compute the force vector (and derivative) for the damping force. */
    Eigen::VectorXd force_damping(const Eigen::VectorXd& config) const;
    Eigen::MatrixXd df_damping(const Eigen::VectorXd& config) const;

    /** Compute the sparse mass matrix used in the kinetic energy computation. */
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mass_matrix() const;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> inv_mass_matrix() const;

    /** Computes the force in the given configuration (according to the enabled forces). */
    Eigen::VectorXd force(const Eigen::VectorXd& config) const;

    /** Computes the derivative of force in the given configuration (according to enabled forces). */
    Eigen::MatrixXd dforce(const Eigen::VectorXd& config) const;

private:
    int32_t particle_unique_id_;
    std::shared_ptr<SimParameters> params_;
    double time_;
    std::vector<Particle, Eigen::aligned_allocator<Particle>> particles_;
    std::vector<Connector*> connectors_;
    std::vector<Saw> saws_;
    std::vector<BendingStencil> bendingStencils_;

    double getTotalParticleMass(int idx) const;

    /**
     * Returns a 2x(2*len(particles)) selection matrix for selecting a given particle.
     * TODO: Eventually optimize this out.
     */
    Eigen::MatrixXd selection_matrix(int index) const;
};

}

#endif
