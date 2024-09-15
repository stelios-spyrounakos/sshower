#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <iostream>

#include "showering.hpp"
#include "kinematic_reco.hpp"

#include <Eigen/Dense>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;
using namespace Eigen;

// returns the unit vector in the direction of the given vector (3D)
Vector3d unit_vector(const Vector3d& vector) {
    return vector / vector.norm();
}

// returns the angle between to vectors (3D)
double angle_between(const Vector3d& v1, const Vector3d& v2) {
    Vector3d v1_unit = unit_vector(v1);
    Vector3d v2_unit = unit_vector(v2);
    double dot_prod = v1_unit.dot(v2_unit);
    if (dot_prod < -1.0) { dot_prod = -1; }
    else if (dot_prod > 1.0) { dot_prod = 1; }
    return acos(dot_prod);
}

// returns the rotation matrix to rotate v1 to the direction of v2
// (uses Rodrigue's Rotation formula)
Matrix3d Get_rotation_matrix(const Vector3d& v1, const Vector3d& v2) {
    Vector3d k = v1.cross(v2);
    Matrix3d I = Matrix3d::Identity();
    if (k.norm() < 1E-12) { return I; }
    Vector3d k_unit = unit_vector(k);
    double angle = angle_between(v1, v2);
    Matrix3d K {
            {0, -k_unit(2), k_unit(1)},
            {k_unit(2), 0, -k_unit(0)},
            {-k_unit(1), k_unit(0), 0}
    };
    // the rotation matrix
    Matrix3d R = I + K * sin(angle) + K * K * (1.0 - cos(angle));
    return R;
}

// rotate 3D vector according to the rotation matrix
Vector3d rotate(const Vector3d& v, const Matrix3d R) {
    Vector3d v_rotated = R * v;
    return v_rotated;
}

// rotates the momenta of the particles (as a whole "jet") from the axis of the
// shower evolution (z-axis) to the direction of the parent in the "lab" frame
vector<Particle> rotate_momenta_lab(const Particle& p,
                                    const vector<Particle>& particles) {
    vector<Particle> rot_particles;
    rot_particles.reserve(particles.size());
    double p_pmag = sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
    Vector3d p_3mom(p.px, p.py, p.pz);
    // the momentum of the parent aligned to the z axis
    Vector3d p_zaxis(0.0, 0.0, p_pmag);
    Matrix3d R = Get_rotation_matrix(p_zaxis, p_3mom);
    Vector3d temp_3mom;
    Vector3d temp_3mom_new;

    vector<Particle>::const_iterator iter;
    for (iter = particles.begin(); iter != particles.end(); ++iter) {
        temp_3mom = {iter->px, iter->py, iter->pz};
        temp_3mom_new = R * temp_3mom;
        rot_particles.push_back(Particle{iter->id, iter->status,
                                            temp_3mom_new(0),
                                            temp_3mom_new(1),
                                            temp_3mom_new(2),
                                            iter->E, iter->m,
                                            iter->t_at_em,
                                            iter->z_at_em,
                                            iter->stopped_evolving});
    }

    return rot_particles;
}

// rotates the momenta of the jet particles with jet 4-momentum q to the
// direction of the progenitor p
vector<Particle> rotate_momenta_prog(const Particle& p,
                                    const vector<double>& q,
                                    const vector<Particle>& jet_particles) {
    vector<Particle> rot_particles;
    rot_particles.reserve(jet_particles.size());
    Vector3d p_3mom(p.px, p.py, p.pz);
    Vector3d q_3mom(q[0], q[1], q[2]);
    Matrix3d R = Get_rotation_matrix(q_3mom, p_3mom);
    Vector3d temp_3mom;
    Vector3d temp_3mom_new;

    vector<Particle>::const_iterator iter;
    for (iter = jet_particles.begin(); iter != jet_particles.end(); ++iter) {
        temp_3mom = {iter->px, iter->py, iter->pz};
        temp_3mom_new = R * temp_3mom;
        rot_particles.push_back(Particle{iter->id, iter->status,
                                            temp_3mom_new(0),
                                            temp_3mom_new(1),
                                            temp_3mom_new(2),
                                            iter->E, iter->m,
                                            iter->t_at_em,
                                            iter->z_at_em,
                                            iter->stopped_evolving});
    }

    return rot_particles;
}

// get the beta factor for the boost for global momentum conservation
// (adapted dorm the Herwig 7 Q-tilde shower). The 3 momenta need to be
// collinear (the rotation must be perfomed before).
// newq is the 4-momentum of the outgoing jet and oldp of the progenitor
vector<double> Get_boost_beta(double k, const vector<double>& newq,
                                const vector<double>& oldp) {
    double q2 = newq[0] * newq[0] + newq[1] * newq[1] + newq[2] * newq[2];
    double q = sqrt(q2);
    double Q2 = newq[3] * newq[3]
                -newq[0] * newq[0] - newq[1] * newq[1] - newq[2] * newq[2];
    double kp = k * sqrt(oldp[0] * oldp[0] + oldp[1] * oldp[1]
                            + oldp[2] * oldp[2]);
    double kp2 = kp * kp;

    // betam = (q*sqrt(qs + Q2) - kp*sqrt(kps + Q2))/(kps + qs + Q2)
    double beta_mag = (q * newq[3]-kp * sqrt(kp2 + Q2)) / (kp2 + q2 + Q2);

    // usually we take the minus sign, since this boost will be smaller.
    // we only require |k \vec p| = |\vec q'| which leaves the sign of
    // the boost open but the 'minus' solution gives a smaller boost
    // parameter, i.e. the result should be closest to the previous
    // result. this is to be changed if we would get many momentum
    // conservation violations at the end of the shower from a hard
    // process.
    Vector3d oldp_3mom(oldp[0], oldp[1], oldp[2]);
    Vector3d beta_temp = -beta_mag * (k / kp) * oldp_3mom;
    vector<double> beta {beta_temp[0], beta_temp[1], beta_temp[2]};
    vector<double> zero {0.0, 0.0, 0.0};
    if (beta_mag >= 1e-12) { return beta; }
    else if (beta_mag >= 1 - 1e-12) {
        beta_mag = 1 - 1e-10;
        beta_temp = -beta_mag * (k / kp) * oldp_3mom;
        beta = {beta_temp[0], beta_temp[1], beta_temp[2]};
        return beta;
    }
    else { return zero; }
}

// apply boost to particles
vector<Particle> boost(const vector<Particle>& particles,
                        const vector<double>& beta) {
    vector<Particle> boosted = particles;
    double betax = beta[0];
    double betay = beta[1];
    double betaz = beta[2];
    double beta_mag = sqrt(betax * betax + betay * betay + betaz * betaz);
    double gamma = 1.0 / sqrt(1.0 - beta_mag * beta_mag);

    int number_of_particles = particles.size();
    for (int i = 0; i < number_of_particles; ++i) {
    double vecx = particles[i].px;
    double vecy = particles[i].py;
    double vecz = particles[i].pz;
    double vec0 = particles[i].E;

    boosted[i].px = - gamma * betax * vec0
            + (1.0 + (gamma - 1.0) * betax * betax / beta_mag / beta_mag) * vecx
            + (gamma - 1.0) * betax * betay * vecy / beta_mag / beta_mag
            + (gamma - 1.0) * betax * betaz * vecz / beta_mag / beta_mag;
    boosted[i].py = - gamma * betay * vec0
            + (gamma - 1.0) * betay * betax * vecx / beta_mag / beta_mag
            + (1.0 + (gamma - 1.0) * betay * betay / beta_mag / beta_mag) * vecy
            + (gamma - 1.0) * betay * betaz * vecz / beta_mag / beta_mag;
    boosted[i].pz = - gamma * betaz * vec0
            + (gamma - 1.0) * betaz * betax * vecx / beta_mag / beta_mag
            + (gamma - 1.0) * betaz * betay * vecy / beta_mag / beta_mag
            + (1.0 + (gamma - 1.0) * betaz * betaz / beta_mag / beta_mag) * vecz;
    boosted[i].E = gamma * (vec0 - betax * vecx - betay * vecy - betaz * vecz);
    }

    return boosted;
}

// the function to solve for obtaining the k rescale factor for mom conservation
// (parameters defined in the header)
double k_fac_eq(double k, void *params) {
    k_fac_eq_params *p = (k_fac_eq_params *)params;
    vector<double> pj2 = p->pj2;
    vector<double> qj2 = p->qj2;
    double sqrt_s_hat = p->sqrt_s_hat;
    double res = 0.0;
    
    int size = pj2.size();
    for (int i = 0; i < size; ++i) {
        res += sqrt(k * k * pj2[i] + qj2[i]);
    }
    res -= sqrt_s_hat;

    return res;

}

// solves for the k factor
// WARNING: problem with the root finding method
double Get_k_fac(vector<double> pj2, vector<double> qj2, double sqrt_s_hat) {
    double tolerance = 0.0001;
    int status;
    double root = 0;
    double k_min = 0.0;
    double k_max = 10.0;
    int iter = 0;
    int max_iter = 1000;
    k_fac_eq_params params = {pj2, qj2, sqrt_s_hat};
    
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    F.function = &k_fac_eq;
    F.params = &params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, k_min, k_max);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        root = gsl_root_fsolver_root(s);
        k_min = gsl_root_fsolver_x_lower(s);
        k_max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(k_min, k_max, 0, tolerance);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return root;
}

// solves for the k factor analytically for 2 jets/progenitors specifically
// (assuming that 3-momentum is properly conserved for the hard process, that
// is, the progenitors)
double Get_k_fac_2_jets(vector<double> pj2, vector<double> qj2,
                    double sqrt_s_hat) {
    if (pj2.size() > 2) {
        cerr << "The function for 2 jets was used when there were actually more"
            << endl;
    }
    double pj2_jet1 = pj2[0];
    double qj2_jet1 = qj2[0];
    double qj2_jet2 = qj2[1];
    double ksq = ( (sqrt_s_hat * sqrt_s_hat - qj2_jet1 - qj2_jet2)
                * (sqrt_s_hat * sqrt_s_hat - qj2_jet1 - qj2_jet2)
                - 4 * qj2_jet1 * qj2_jet2 )
                / (4 * sqrt_s_hat * sqrt_s_hat *pj2_jet1);
    if (ksq < 1e-12) { ksq = 1e-12;}
    double k = sqrt(ksq);
    return k;
}

// implements global momentum conservation
vector<Particle> global_mom_cons(const vector<Jet>& showered_jets) {
    // initializing total energy
    double sqrt_s_hat = 0.0;
    // pj are the progenitor 3-momenta and qj are the showered jet 4-momenta
    vector<double> pj2;
    vector<double> qj2;
    // here oldp is the progenitor momenta and newq the showered jet momenta as
    // in the naming scheme for the boost
    vector< vector<double> > oldps;
    vector< vector<double> > newqs;
    vector<double> oldp {0.0, 0.0, 0.0, 0.0};
    vector<double> newq {0.0, 0.0, 0.0, 0.0};
    vector<double> q_temp {0.0, 0.0, 0.0, 0.0};

    // iterate through the jets
    int number_of_jets = showered_jets.size();

    for (int i = 0; i < number_of_jets; ++i) {
        // reset newq
        newq = {0.0, 0.0, 0.0, 0.0};
        Jet jet = showered_jets[i];
        oldp = {jet.progenitor.px, jet.progenitor.py,
                jet.progenitor.pz, jet.progenitor.E};
        pj2.push_back(oldp[0] * oldp[0] + oldp[1] * oldp[1]
                        + oldp[2] * oldp[2]);
        oldps.push_back(oldp);
        sqrt_s_hat += oldp[3];

        // iterate through the particles of the jet
        int number_of_jetparts = jet.particles.size();
        for (int j = 0; j < number_of_jetparts; ++j) {
            q_temp = {jet.particles[j].px, jet.particles[j].py,
                        jet.particles[j].pz, jet.particles[j].E};
            newq[0] += q_temp[0];
            newq[1] += q_temp[1];
            newq[2] += q_temp[2];
            newq[3] += q_temp[3];
        }

        double qj2_temp = newq[3] * newq[3] - newq[0] * newq[0]
                            - newq[1] * newq[1] - newq[2] * newq[2];
        if (isnan(qj2_temp)) { qj2_temp = 0.0; }
        qj2.push_back(qj2_temp);
        newqs.push_back(newq);
    }

    ///// FOR DEBUGGING /////
    for (int i = 0; i < number_of_jets; ++i) {
        vector<double> oldpp = oldps[i];
        vector<double> newqq = newqs[i];
        cout << "JET " << i << " PROGENITOR pj2: " << pj2[i]
            << ", JET qj2: " << qj2[i] << endl;

        cout << "JET " << i << " PROGENITOR oldp: (" << oldpp[0] << ", "
            << oldpp[1] << ", " << oldpp[2] << ", " << oldpp[3] << ")" << endl;

        cout << "JET " << i << " JET newq: (" << newqq[0] << ", "
            << newqq[1] << ", " << newqq[2] << ", " << newqq[3] << ")" << endl;
    }
    /////////////////////////

    // get the k factor for the boosts
    double k = Get_k_fac_2_jets(pj2, qj2, sqrt_s_hat);

    ///// FOR DEBUGGING /////
    cout << "K FACTOR: " << k << endl;
    /////////////////////////

    // now we rotate and boost the particles in the jets for global momentum
    // consrvation
    vector<Particle> showered_particles_boosted;
    int total_jetparticles = 0;

    ///// FOR DEBUGGING /////
    cout << "JET RECO LOOPS START" << k << endl;
    /////////////////////////

    for (int i = 0; i < number_of_jets; ++i) {
       total_jetparticles += showered_jets.at(i).particles.size();    
    }
    showered_particles_boosted.reserve(total_jetparticles);

    for (int i = 0; i < number_of_jets; ++i) {
        Jet jet = showered_jets.at(i);
        
        vector<Particle> rotated;
        vector<double> beta;
        vector<Particle> boosted;

        int number_of_jetparts = jet.particles.size();
        if (number_of_jetparts == 1) {
            // do nothing if the particle has not showered
            showered_particles_boosted.push_back(jet.particles.at(0));
        }
        else {
            // rotate and boost the jet particles
            rotated = rotate_momenta_prog(jet.progenitor, newqs.at(i),
                                        jet.particles);

            ///// FOR DEBUGGING /////
            cout << "ROTATION SUCCESS" << endl;
            /////////////////////////

            ///// FOR DEBUGGING /////
            cout << "COMPUTING BOOST BETA FACTOR: " << endl;
            /////////////////////////

            beta = Get_boost_beta(k, newqs.at(i), oldps.at(i));
            
            ///// FOR DEBUGGING /////
            cout << "BOOST BETA FACTOR: "
                << beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2]
                << endl;
            /////////////////////////

            boosted = boost(rotated, beta);

            ///// FOR DEBUGGING /////
            cout << "BOOST SUCCESS" << endl;
            /////////////////////////

            ///// FOR DEBUGGING /////
            cout << "JET PARTICLES RECO LOOP START" << k << endl;
            /////////////////////////
            for (int j = 0; j < number_of_jetparts; ++j) {
                showered_particles_boosted.push_back(boosted.at(j));
            }
            ///// FOR DEBUGGING /////
            cout << "JET PARTICLES RECO LOOP END" << k << endl;
            /////////////////////////

            ///// FOR DEBUGGING /////
            cout << "RECO PARTICLES ADDED TO VECTOR" << endl;
            /////////////////////////
        }
    }

    ///// FOR DEBUGGING /////
    cout << "JET RECO LOOPS END" << k << endl;
    /////////////////////////

    return showered_particles_boosted;
}

// returns the total 3-momentum of an event
std::vector<double> check_mom_cons(Event event) {
    vector<double> total(3, 0.0);
    for (int i = 0; i < event.particles.size(); ++i) {
        Particle particle = event.particles.at(i);
        if (particle.status == 1) {
            total[0] += particle.px;
            total[1] += particle.py;
            total[2] += particle.pz;
        }
    }
    return total;
}