#ifndef GUARD_shower_helpers_hpp
#define GUARD_shower_helpers_hpp

#include "splittings.hpp"
#include "alphaS.hpp"

// returns alphaS at the emission scale (pT of the emission)
double Get_alphaS(double t, double z, double Q_cutoff, const alphaS& aS);

// returns the overestimated value of alphaS (at the cutoff scale)
double Get_alphaS_over(double Q_cutoff, const alphaS& aS);

// they return the (overestimated) upper and lower z limits
double zp_over(double t, double Q_cutoff);

double zm_over(double t, double Q_cutoff);

// returns the momentum fraction candidate of the emission
// for the appropriate splitting
double Get_z_em(double (*inverse_rho)(double, double),
                double (*rho)(double, double),
                double t, double Q_cutoff, double aS_over, double Rand);

// returns the transverse momentum squared of the emission
double Get_pT2(double t, double z);

// returns the virtual mass squared of the emitter
double Get_mvirt2(double t, double z);

// the function to solve for obtaining the evolution scale of the emission
double em_scale_func(double t, void *params);

// the sctruct for the parameters for the function
struct em_scale_func_params
{
    double Q;
    double Q_cutoff;
    double aS_over;
    double Rand;
    int branching_type;    // 1 for q->qg, 2 for g->qqbar, 3 for g->gg
};

// returns the emission scale and modifies the bool continue_evol
double Get_t_em(double Q, double Q_cutoff, double aS_over, double Rand,
                double t_fac, int branching_type, bool& continue_evol);

#endif