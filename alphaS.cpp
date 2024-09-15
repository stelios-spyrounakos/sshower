#define _USE_MATH_DEFINES

#include <cmath>
#include "constants.hpp"
#include "alphaS.hpp"

using namespace std;

// constructor
alphaS::alphaS(double asmz, double mz, double mb, double mc, int order)
    : order(order), asmz(asmz), mz(mz), mz2(mz * mz), mb2(mb * mb), mc2(mc * mc) {
    asmb = alphasQ(mb);
    asmc = alphasQ(mc);
}

// function to access alphaS at scale Q
double alphaS::alphasQ(double Q) const {
    if (order == 0) {
        return As0(Q * Q);
    } else {
        return As1(Q * Q);
    }
}

// NLO alphaS:
double alphaS::As1(double t) const {
    double tref, asref, b0, b1;
    if (t >= mb2) {
        tref = mz2;
        asref = asmz;
        b0 = Beta0(5) / (2.0 * M_PI);
        b1 = Beta1(5) / (2.0 * M_PI * 2.0 * M_PI);
    } else if (t >= mc2) {
        tref = mb2;
        asref = asmb;
        b0 = Beta0(4) / (2.0 * M_PI);
        b1 = Beta1(4) / (2.0 * M_PI * 2.0 * M_PI);
    } else {
        tref = mc2;
        asref = asmc;
        b0 = Beta0(3) / (2.0 * M_PI);
        b1 = Beta1(3) / (2.0 * M_PI * 2.0 * M_PI);
    }
    double w = 1.0 + b0 * asref * log(t / tref);
    return asref / w * (1.0 - b1 / b0 * asref * log(w) / w);
}

// LO alphaS:
double alphaS::As0(double t) const {
    double tref, asref, b0;
    if (t >= mb2) {
        tref = mz2;
        asref = asmz;
        b0 = Beta0(5) / (2.0 * M_PI);
    } else if (t >= mc2) {
        tref = mb2;
        asref = asmb;
        b0 = Beta0(4) / (2.0 * M_PI);
    } else {
        tref = mc2;
        asref = asmc;
        b0 = Beta0(3) / (2.0 * M_PI);
    }
    return 1.0 / (1.0 / asref + b0 * log(t / tref));
}

// beta functions:
double alphaS::Beta0(int nf) const {
    return (11.0 / 6.0 * CA) - (2.0 / 3.0 * TR * nf);
}

double alphaS::Beta1(int nf) const {
    return (17.0 / 6.0 * CA * CA) - ((5.0 / 3.0 * CA + CF) * TR * nf);
}