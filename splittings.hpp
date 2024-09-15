#ifndef GUARD_splittings_hpp
#define GUARD_splittings_hpp

// q -> qg splitting function
double Pqq(double z);

// q -> qg splitting function overestimate
double Pqq_over(double z);

// q -> qg rho tilde function
double rho_qq(double z, double aS_over);

// q -> qg rho tilde function's inverse
double inverse_rho_qq(double z, double aS_over);



// g -> qqbar splitting function
double Pgq(double z);

// g -> qqbar splitting function overestimate
double Pgq_over(double z);

// g -> qqbar rho tilde function
double rho_gq(double z, double aS_over);

// g -> qqbar rho tilde function's inverse
double inverse_rho_gq(double z, double aS_over);



// g -> gg splitting function
double Pgg(double z);

// g -> gg splitting function overestimate
double Pgg_over(double z);

// g -> gg rho tilde function
double rho_gg(double z, double aS_over);

// g -> gg rho tilde function's inverse
double inverse_rho_gg(double z, double aS_over);

#endif