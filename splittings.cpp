#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>

#include "constants.hpp"
#include "splittings.hpp"

using namespace std;

// q -> qg splitting function
double Pqq(double z) {
    return CF * (1.0 + z * z) / (1.0 - z);
}

// q -> qg splitting function overestimate
double Pqq_over(double z) {
    return CF * 2 / (1.0 - z);
}

// q -> qg rho tilde function
double rho_qq(double z, double aS_over) {
    return -2.0 * CF * (aS_over / (2.0 * M_PI)) * log1p(-z);
}

// q -> qg rho tilde function's inverse
double inverse_rho_qq(double z, double aS_over) {
    return 1 - exp(-0.5 * z / (CF * aS_over / (2.0 * M_PI)));
}



// g -> qqbar splitting function
double Pgq(double z) {
    return TR * (1.0 - 2.0 * z * (1.0 - z));
}

// g -> qqbar splitting function overestimate
double Pgq_over(double z) {
    return TR;
}

// g -> qqbar rho tilde function
double rho_gq(double z, double aS_over) {
    return TR * (aS_over / (2.0 * M_PI)) * z;
}

// g -> qqbar rho tilde function's inverse
double inverse_rho_gq(double z, double aS_over) {
    return (1.0 / (TR * aS_over / (2.0 * M_PI))) * z;
}



// g -> gg splitting function
double Pgg(double z) {
    return CA * (z / (1.0 - z) + (1.0 - z) / z + z * (1.0 - z));
}

// g -> gg splitting function overestimate
double Pgg_over(double z) {
    return CA * (1.0 / (1.0 - z) + 1.0 / z);
}

// g -> gg rho tilde function
double rho_gg(double z, double aS_over) {
    return -CA * (aS_over / (2.0 * M_PI)) * log(1.0 / z - 1.0);
}

// g -> gg rho tilde function's inverse
double inverse_rho_gg(double z, double aS_over) {
    return 1 / (1 + exp(- z / (CA * aS_over / (2.0 * M_PI))));
}