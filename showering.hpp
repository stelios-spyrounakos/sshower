#ifndef GUARD_showering_hpp
#define GUARD_showering_hpp

#include <vector>

// struct that contains emission information
struct Emission_info {
    double t_em = 0;
    double z_em = 1;
    double pT2_em = 0;
    double phi = 0;
    double mvirt2 = 0;
    bool continue_evolution = false;
    bool generated_em = false;
};

// generates emission information
Emission_info generate_emission(double Q, double Q_cutoff, double t_fac,
                                unsigned int branching_type);

struct Particle {
    int id;
    int status;
    double px;
    double py;
    double pz;
    double E;
    double m;
    // store the t and z when the particle is emitted. To be used for the
    // scale of further emissions
    double t_at_em;
    double z_at_em;
    bool stopped_evolving;
};

struct Jet {
    Particle progenitor;
    std::vector<Particle> particles;
};

struct Event {
    std::vector<Particle> particles;
};

// perform a single 1->2 branching (pb gets the z fraction and pc the 1-z)
void branching_1to2(Particle& pa, Particle& pb, Particle& pc,
                    double Q_cutoff, double t_fac, double t_cutoff);

// perform a single 1->2 branching (pb gets the z fraction and pc the 1-z)
void branching_1to2_v2(Particle& pa, Particle& pb, Particle& pc,
                    double Q_cutoff, double t_fac, double t_cutoff);

// shower a particle
std::vector<Particle> shower_particle(const Particle& p, double Q_cutoff,
                                    double t_fac, double t_cutoff);

// shower a particle (p is in the emissions vector)
void shower_particle_v2(Particle& p, double Q_cutoff, double t_fac,
                    double t_cutoff, std::vector<Particle>& emissions);

// shower a progenitor
Jet shower_progenitor(const Particle& prog, double Q_cutoff, double t_fac,
                    double t_cutoff);

// shower an event (showers the appropriate hard process generated particles)
std::vector<Particle> shower_event(Event& event, double Q_cutoff);

#endif