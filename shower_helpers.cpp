#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "splittings.hpp"
#include "alphaS.hpp"
#include "shower_helpers.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;

// returns alphaS at the emission scale (pT of the emission)
double Get_alphaS(double t, double z, double Q_cutoff, const alphaS& aS) {
    double scale = z * (1.0 - z) * sqrt(t);     // = pT_emission
    if (scale < Q_cutoff) { return aS.alphasQ(Q_cutoff); }
    return aS.alphasQ(scale);
}

// returns the overestimated value of alphaS (at the cutoff scale)
double Get_alphaS_over(double Q_cutoff, const alphaS& aS) {
    return aS.alphasQ(Q_cutoff);
}

// they return the (overestimated) upper and lower z limits
double zp_over(double t, double Q_cutoff) {
    return 1.0 - sqrt(Q_cutoff * Q_cutoff / t);
}

double zm_over(double t, double Q_cutoff) {
    return sqrt(Q_cutoff * Q_cutoff / t);
}

// returns the momentum fraction candidate of the emission
// for the appropriate splitting
double Get_z_em(double (*inverse_rho)(double, double),
                double (*rho)(double, double),
                double t, double Q_cutoff, double aS_over, double Rand) {
    double arg = rho(zm_over(t, Q_cutoff), aS_over) + Rand * 
                    (rho(zp_over(t, Q_cutoff), aS_over) -
                    rho(zm_over(t, Q_cutoff), aS_over));
    return inverse_rho(arg, aS_over);
}

// returns the transverse momentum squared of the emission
double Get_pT2(double t, double z) {
    return z * z * (1.0 - z) * (1.0 - z) * t;
}

// returns the virtual mass squared of the emitter
double Get_mvirt2(double t, double z) {
    return z * (1.0 - z) * t;
}

// the function to solve for obtaining the evolution scale of the emission
// (parameters defined in the header)
double em_scale_func(double t, void *params) {
    em_scale_func_params *p = (em_scale_func_params *)params;
    double Q = p->Q;
    double Q_cutoff = p->Q_cutoff;
    double aS_over = p->aS_over;
    double Rand = p->Rand;
    int branching_type = p->branching_type;

    double r = 0.0;
    
    // 1 for q->qg, 2 for g->qqbar, 3 for g->gg
    switch (branching_type) {
    case 1:
        r = rho_qq(zp_over(t, Q_cutoff), aS_over)
                    - rho_qq(zm_over(t, Q_cutoff), aS_over);
        break;
    case 2:
        r = rho_gq(zp_over(t, Q_cutoff), aS_over)
                    - rho_gq(zm_over(t, Q_cutoff), aS_over);
        break;
    case 3:
        r = rho_gg(zp_over(t, Q_cutoff), aS_over)
                    - rho_gg(zm_over(t, Q_cutoff), aS_over);
        break;
    default:
        cerr << "Invalid branching" << endl;
        break;
    }

    return log(t / Q / Q) - (1.0 / r) * log(Rand);
}

// returns the emission scale and modifies the bool continue_evol
// (uses the GSL's implementation of the Brent-Dekker algorithm)
// t_fac is a factor that modifies the allowed range of t, so a t too close to
// the cutoff t_c is not returned (t_c = 4*t_0, where t_0 = Q_cutoff^2)
double Get_t_em(double Q, double Q_cutoff, double aS_over, double Rand,
                double t_fac, int branching_type, bool& continue_evol) {
    double tolerance = 0.0001;
    int status;
    double root = 0;
    double t_min = t_fac * Q_cutoff * Q_cutoff;
    double t_max = Q * Q;
    
    ///// FOR DEBUGGING /////
    cout << "FINDING T_EM BETWEEN " << t_min << " AND " << t_max << endl;
    /////////////////////////

    int iter = 0;
    int max_iter = 1000;
    em_scale_func_params params = {Q, Q_cutoff, aS_over, Rand, branching_type};
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    F.function = &em_scale_func;
    F.params = &params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, t_min, t_max);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        root = gsl_root_fsolver_root(s);
        t_min = gsl_root_fsolver_x_lower(s);
        t_max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(t_min, t_max, 0, tolerance);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);

    if ( abs(em_scale_func(root, &params)) > tolerance ) {
        continue_evol = false;
    }
    else if (root > Q * Q) {
        continue_evol = false;
    }
    else if (root < t_fac * Q_cutoff * Q_cutoff) {
        continue_evol = false;
    }
    else {
        continue_evol = true;
    }

    return root;
}