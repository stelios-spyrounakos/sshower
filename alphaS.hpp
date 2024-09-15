#ifndef GUARD_alphaS_hpp
#define GUARD_alphaS_hpp

class alphaS {
public:
    // constructor
    alphaS(double asmz, double mz, double mb = 4.75, double mc = 1.27, int order = 1);

    // function to access alphaS at scale Q
    double alphasQ(double Q) const;

private:
    // data members
    int order;
    double asmz;
    double mz;
    double mz2;
    double mb2;
    double mc2;
    double asmb;
    double asmc;

    // NLO alphaS:
    double As1(double t) const;

    // LO alphaS:
    double As0(double t) const;

    // beta functions:
    double Beta0(int nf) const;
    double Beta1(int nf) const;
};

#endif