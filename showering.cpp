#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
//#include <iterator>
//#include <ctime>
#include "splittings.hpp"
#include "alphaS.hpp"
#include "shower_helpers.hpp"
#include "showering.hpp"
#include "kinematic_reco.hpp"

using namespace std;

// generates emission information
Emission_info generate_emission(double Q, double Q_cutoff, double t_fac,
                                int branching_type) {
    // set up the random number generator
    //static mt19937 gen(time(nullptr));
    static random_device rd;        // to obtain a random seed
    static mt19937 gen(rd());       // standard mersenne twister engine
    static uniform_real_distribution<> unif(0.0000000001, 1.0);
    // produce the random numbers
    double Rand1 = unif(gen);
    double Rand2 = unif(gen);
    double Rand3 = unif(gen);
    double Rand4 = unif(gen);

    // initialize alphaS object
    static alphaS aS(0.118, 91.1876);
    static double aS_over = Get_alphaS_over(Q_cutoff, aS);
    // proceed with generating the emission info
    Emission_info em_info;
    em_info.generated_em = true;
    // get the (candidate) emission scale

    ///// FOR DEBUGGING /////
    cout << "GENERATING t_em FOR BRANCHING TYPE " << branching_type << endl;
    /////////////////////////

    bool cont_evol;
    em_info.t_em = Get_t_em(Q, Q_cutoff, aS_over, Rand1, t_fac,
                    branching_type, cont_evol);
    em_info.continue_evolution = cont_evol;

    ///// FOR DEBUGGING /////
    cout << "GENERATED t_em: " << em_info.t_em << " AND continue_evolution: "
        << em_info.continue_evolution << endl;
    /////////////////////////

    // stop if no solution is found
    if (em_info.continue_evolution == false) {
        em_info.t_em = 0.0;
        em_info.z_em = 1.0;
        em_info.pT2_em = 0.0;
        em_info.phi = 0.0;
        em_info.mvirt2 = 0.0;
        em_info.generated_em = false;
        return em_info;
    }

    // get the (candidate) mom fraction of the emission for the aproppriate
    // branching type (1 for q->qg, 2 for g->qqbar, 3 for g->gg)
    switch (branching_type) {
    case 1:
        em_info.z_em = Get_z_em(inverse_rho_qq, rho_qq, em_info.t_em, Q_cutoff,
                                aS_over, Rand2);
        break;
    case 2:
        em_info.z_em = Get_z_em(inverse_rho_gq, rho_gq, em_info.t_em, Q_cutoff,
                                aS_over, Rand2);
        break;
    case 3:
        em_info.z_em = Get_z_em(inverse_rho_gg, rho_gg, em_info.t_em, Q_cutoff,
                                aS_over, Rand2);
        break;
    default:
        cerr << "Invalid branching" << endl;
        break;
    }

    // get the transverse momentum squared of the (candidate) emission
    em_info.pT2_em = Get_pT2(em_info.t_em, em_info.z_em);
    // generate angle for the pT of the emission
    // (uniformly, does not take into account spin correlations)
    em_info.phi = (2 * unif(gen) - 1) * M_PI;

    ///// FOR DEBUGGING /////
    cout << "VETOES START" << endl;
    /////////////////////////

    // perform the vetos
    if (em_info.pT2_em < 0.0) {
        em_info.generated_em = false;
        ///// FOR DEBUGGING /////
        cout << "VETOED DUE TO pT2" << endl;
        /////////////////////////
    }
    switch (branching_type) {
    case 1:
        if (Pqq(em_info.z_em) / Pqq_over(em_info.z_em) < Rand3) {
            em_info.generated_em = false;
            ///// FOR DEBUGGING /////
            cout << "VETOED DUE TO KERNEL" << endl;
            /////////////////////////
        }
        break;
    case 2:
        if (Pgq(em_info.z_em) / Pgq_over(em_info.z_em) < Rand3) {
            em_info.generated_em = false;
            ///// FOR DEBUGGING /////
            cout << "VETOED DUE TO KERNEL" << endl;
            /////////////////////////
        }
        break;
    case 3:
        if (Pgg(em_info.z_em) / Pgg_over(em_info.z_em) < Rand3) {
            em_info.generated_em = false;
            ///// FOR DEBUGGING /////
            cout << "VETOED DUE TO KERNEL" << endl;
            /////////////////////////
        }
        break;
    default:
        cerr << "Invalid branching" << endl;
        break;
    }
    if (Get_alphaS(em_info.t_em, em_info.z_em, Q_cutoff, aS)
        / aS_over < Rand4) {
        em_info.generated_em = false;
        ///// FOR DEBUGGING /////
        cout << "VETOED DUE TO aS" << endl;
        /////////////////////////
    }

    ///// FOR DEBUGGING /////
    cout << "VETOES OVER" << endl;
    /////////////////////////

    // get the virtual mass squared
    em_info.mvirt2 = Get_mvirt2(em_info.t_em, em_info.z_em);

    if (em_info.generated_em == false) {
        em_info.z_em = 1.0;
        em_info.pT2_em = 0.0;
        em_info.phi = 0.0;
        em_info.mvirt2 = 0.0;
    }
    
    return em_info;
}

// perform a single 1->2 branching (pb gets the z fraction and pc the 1-z)
void branching_1to2(Particle& pa, Particle& pb, Particle& pc,
                    double Q_cutoff, double t_fac, double t_cutoff) {
    // find the emission type and generate emission information
    int branching_type;
    Emission_info em_info;
    Emission_info em_info_temp;
    double Q = pa.z_at_em * sqrt(pa.t_at_em);
    
    // try branching while the emission is vetoed. Stop if an emission is
    // generated or the evolution of the particle cannot continue
    do {
        if (Q < sqrt(t_cutoff)) { 
            pa.stopped_evolving = true;
            return; 
        }
        // for quarks we only have q->qg emissions (we ignore the top) and the
        // flavor of the "produced" quark is the same
        if (abs(pa.id) > 0 && abs(pa.id) < 6) {
            branching_type = 1;
            em_info = generate_emission(Q, Q_cutoff, t_fac, branching_type);
            // the quark and the gluon
            pb.id = pa.id;
            pc.id = 21;
        }
        // for gluons we pick the emission with the largest t out of the
        // possible ones
        else if (abs(pa.id) == 21) {
            // check for branching type 3 (g->gg) and the 5 (ignoring top)
            // flavor possibilities of branching type 2 (g->qqbar). Start with 3
            branching_type = 3;
            em_info = generate_emission(Q, Q_cutoff, t_fac, branching_type);
            // set t to 0 if no solution is found here so that any solution
            // found next will win
            if (em_info.continue_evolution == false) {
                em_info.t_em = 0;
            }
            // the gluons
            pb.id = 21;
            pc.id = 21;
            // if an emission has a larger t we will replace the previous
            branching_type = 2;
            for (int flavor = 1; flavor < 6; ++flavor) {
                em_info_temp = generate_emission(Q, Q_cutoff, t_fac,
                                                branching_type);
                // skip if no solution for t is found
                if (em_info_temp.continue_evolution == false) {
                    continue;
                }

                if (em_info_temp.t_em > em_info.t_em) {
                    em_info = em_info_temp;
                    // the quark antiquark pair
                    pb.id = flavor;
                    pc.id = -flavor;
                }
            }
        }
        Q = sqrt(em_info.t_em);
    } while (em_info.generated_em == false
            && em_info.continue_evolution == true);

    ///// FOR DEBUGGING /////
    cout << "EM INFO FROM branching_1to2:" << endl;
    cout << "t_em: " << em_info.t_em << ", z_em: " << em_info.z_em
        << ", pT2_em: " << em_info.pT2_em << ", phi: " << em_info.phi
        << ", mvirt2: " << em_info.mvirt2 
        << ", continue_evolution: " << em_info.continue_evolution
        << ", generated_em: " << em_info.generated_em << endl;
    /////////////////////////

    if (em_info.continue_evolution == false) {
        pa.stopped_evolving = true;
        return;
    }
    // t and z for b and c at emission
    pb.t_at_em = em_info.t_em;
    pc.t_at_em = em_info.t_em;
    pb.z_at_em = em_info.z_em;
    pc.z_at_em = 1.0 - em_info.z_em;

    // constructing the momenta of pb and pc for pa in the z axis
    double pT = sqrt(em_info.pT2_em);
    double pmag = sqrt(pa.px * pa.px + pa.py * pa.py + pa.pz * pa.pz);

    pb.px = pT * cos(em_info.phi);
    pb.py = pT * sin(em_info.phi);
    pb.pz = em_info.z_em * pmag;
    pb.E = sqrt(pb.px * pb.px + pb.py * pb.py + pb.pz * pb.pz);

    pc.px = - pT * cos(em_info.phi);
    pc.py = - pT * sin(em_info.phi);
    pc.pz = (1.0 - em_info.z_em) * pmag;
    pc.E = sqrt(pc.px * pc.px + pc.py * pc.py + pc.pz * pc.pz);

    // the rest of the particle info
    pb.status = 1;
    pc.status = 1;
    pb.m = 0;
    pc.m = 0;
    pb.stopped_evolving = false;
    pc.stopped_evolving = false;

    // rotate particles b and c to the original direction of their parent a
    vector<Particle> emissions{pb, pc};
    
    ///// FOR DEBUGGING /////
    cout << "PARENT PARTICLE A INFO: " << endl;
    cout << "id: " << pa.id << ", status: " << pa.status << ", px: " << pa.px
        << ", py: " << pa.py << ", pz: " << pa.pz << ", E: " << pa.E
        << ", m: " << pa.m << ", t_at_em: " << pa.t_at_em
        << ", z_at_em: " << pa.z_at_em << ", stopped_evolving: "
        << pa.stopped_evolving << endl;

    cout << "EMITTED PARTICLE B INFO BEFORE ROT: " << endl;
    cout << "id: " << pb.id << ", status: " << pb.status << ", px: " << pb.px
        << ", py: " << pb.py << ", pz: " << pb.pz << ", E: " << pb.E
        << ", m: " << pb.m << ", t_at_em: " << pb.t_at_em
        << ", z_at_em: " << pb.z_at_em << ", stopped_evolving: "
        << pb.stopped_evolving << endl;

    cout << "EMITTED PARTICLE C INFO BEFORE ROT: " << endl;
    cout << "id: " << pc.id << ", status: " << pc.status << ", px: " << pc.px
        << ", py: " << pc.py << ", pz: " << pc.pz << ", E: " << pc.E
        << ", m: " << pc.m << ", t_at_em: " << pc.t_at_em
        << ", z_at_em: " << pc.z_at_em << ", stopped_evolving: "
        << pc.stopped_evolving << endl;
    /////////////////////////

    vector<Particle> emissions_rot = rotate_momenta_lab(pa, emissions);
    pb = emissions_rot[0];
    pc = emissions_rot[1];

    ///// FOR DEBUGGING /////
    cout << "EMITTED PARTICLE B INFO AFTER ROT: " << endl;
    cout << "id: " << pb.id << ", status: " << pb.status << ", px: " << pb.px
        << ", py: " << pb.py << ", pz: " << pb.pz << ", E: " << pb.E
        << ", m: " << pb.m << ", t_at_em: " << pb.t_at_em
        << ", z_at_em: " << pb.z_at_em << ", stopped_evolving: "
        << pb.stopped_evolving << endl;

    cout << "EMITTED PARTICLE C INFO AFTER ROT: " << endl;
    cout << "id: " << pc.id << ", status: " << pc.status << ", px: " << pc.px
        << ", py: " << pc.py << ", pz: " << pc.pz << ", E: " << pc.E
        << ", m: " << pc.m << ", t_at_em: " << pc.t_at_em
        << ", z_at_em: " << pc.z_at_em << ", stopped_evolving: "
        << pc.stopped_evolving << endl;
    /////////////////////////

    return;
}

// shower a progenitor
Jet shower_progenitor(const Particle& p, double Q_cutoff, double t_fac,
                    double t_cutoff) {
    Jet jet;
    jet.progenitor = p;
    jet.particles.reserve(50);
    jet.particles.push_back(p);
    int i = 0;
    // this loops branches the b particle of the previous branch (replaces pa
    // with pb and adds pc to the jet but leaves the counter the same so the
    // new pb becomes pa in the next loop) and adds pc to the jet. When the
    // "chain" of b particles cannot branch anymore, the counter is increased
    // and the process repeats with the first pc added to the jet as pa and
    // when its own now chain of b particles ends, it moves on to the next pc
    // added and the cycle continues

    ///// FOR DEBUGGING /////
    cout << "STARTING SHOWERING PROGENITOR LOOP" << endl;
    /////////////////////////
    while (i < jet.particles.size()) {
        Particle pa = jet.particles.at(i);
        Particle pb = {0, 0, 0, 0, 0, 0, 0, 0, 0, true};
        Particle pc = {0, 0, 0, 0, 0, 0, 0, 0, 0, true};
        branching_1to2(pa, pb, pc, Q_cutoff, t_fac, t_cutoff);
        if (pa.stopped_evolving == true) {
            ++i;
            continue;
        }
        jet.particles.at(i) = pb;
        jet.particles.push_back(pc);
    }
    jet.particles.shrink_to_fit();
    ///// FOR DEBUGGING /////
    cout << "SHOWERING PROGENITOR LOOP END" << endl;
    /////////////////////////

    return jet;
}



// shower an event (showers the appropriate hard process generated particles)
vector<Particle> shower_event(Event& event, double Q_cutoff) {
    // the minimun scale
    double t_min = Q_cutoff * Q_cutoff;
    // the cutoff on the emission scale
    double t_cutoff = 4 * t_min;
    // t_fac is a factor that modifies the allowed range of t, so a t too close
    // to the cutoff t_c is not returned (t_c = 4*t_0, where t_0 = Q_cutoff^2)
    double t_fac = 3.999;

    vector<Particle> final_particles;
    vector<Jet> jets;
    jets.reserve(2);

    // produce the jets
    vector<Particle>::const_iterator p;
    int i = 0;
    for (p = event.particles.begin(); p != event.particles.end(); ++p) {
        if (abs(p->id) == 11) {
            // just adds initial state particles
            final_particles.push_back(*p);
            continue;
        }

        ///// FOR DEBUGGING /////
        cout << "SHOWERING PROGENITOR " << i << endl;
        /////////////////////////

        Jet jet = shower_progenitor(*p, Q_cutoff, t_fac, t_cutoff);

        ///// FOR DEBUGGING /////
        cout << "SHOWERED PROGENITOR " << i << endl;
        /////////////////////////

        i++;
        jets.push_back(jet);
    }
    jets.shrink_to_fit();

    // apply global momentum conservation

    ///// FOR DEBUGGING /////
    cout << "APPLY MOM CON" << endl;
    /////////////////////////

    vector<Particle> mom_con_jet_parts = global_mom_cons(jets);

    ///// FOR DEBUGGING /////
    cout << "MOM CON DONE" << endl;
    ///// FOR DEBUGGING /////

    // merge the initial particles with the jets
    final_particles.reserve(final_particles.size() + mom_con_jet_parts.size());
    
    final_particles.insert(final_particles.end(),
                        make_move_iterator(mom_con_jet_parts.begin()),
                        make_move_iterator(mom_con_jet_parts.end()));
    

    move(mom_con_jet_parts.begin(), mom_con_jet_parts.end(),
        back_inserter(final_particles));

    ///// FOR DEBUGGING /////
    cout << "FINAL PARTICLES ADDED TO VECTOR" << endl;
    /////////////////////////

    return final_particles;
}


