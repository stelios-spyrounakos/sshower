#ifndef GUARD_kinematic_reco_hpp
#define GUARD_kinematic_reco_hpp

#include <vector>
#include "showering.hpp"

#include <Eigen/Dense>

// returns the unit vector in the direction of the given vector (3D)
Eigen::Vector3d unit_vector(const Eigen::Vector3d& vector);

// returns the angle between to vectors (3D)
double angle_between(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

// returns the rotation matrix to rotate v1 to the direction of v2
// (uses Rodrigue's Rotation formula)
Eigen::Matrix3d Get_rotation_matrix(const Eigen::Vector3d& v1,
                                    const Eigen::Vector3d& v2);

// rotate 3D vector according to the rotation matrix
Eigen::Vector3d rotate(const Eigen::Vector3d& v, const Eigen::Matrix3d R);

// rotates the momenta of the particles (as a whole "jet" from the axis of the
// shower evolution (z-axis) to the direction of the parent in the "lab" frame
std::vector<Particle> rotate_momenta_lab(const Particle& p,
                                    const std::vector<Particle>& particles);

// rotates the momenta of the jet particles with jet 4-momentum q to the
// direction of the progenitor p
std::vector<Particle> rotate_momenta_prog(const Particle& p,
                                    const std::vector<double>& q,
                                    const std::vector<Particle>& jet_particles);

// get the beta factor for the boost for global momentum conservation
// (adapted dorm the Herwig 7 Q-tilde shower). The 3 momenta need to be
// collinear (the rotation must be perfomed before).
// newq is the 4-momentum of the outgoing jet and oldp of the progenitor
std::vector<double> Get_boost_beta(double k, const std::vector<double>& newq,
                                    const std::vector<double>& oldp);

// apply boost to particles
std::vector<Particle> boost(const std::vector<Particle>& particles,
                            const std::vector<double>& beta);

// the function to solve for obtaining the k rescale factor for mom conservation
double k_fac_eq(double k, void *params);

// the sctruct for the parameters for the function
struct k_fac_eq_params
{
    std::vector<double> pj2;
    std::vector<double> qj2;
    double sqrt_s_hat;
};

// solves for the k factor
// WARNING: problem with the root finding method
double Get_k_fac(std::vector<double> pj2, std::vector<double> qj2,
                    double sqrt_s_hat);

// solves for the k factor analytically for 2 jets/progenitors specifically
// (assuming that 3-momentum is properly conserved for the hard process, that
// is, the progenitors)
double Get_k_fac_2_jets(std::vector<double> pj2, std::vector<double> qj2,
                    double sqrt_s_hat);

// implements global momentum conservation
std::vector<Particle> global_mom_cons(const std::vector<Jet>& showered_jets);

// returns the total 3-momentum of an event
std::vector<double> check_mom_cons(Event event);

#endif