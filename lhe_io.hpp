#ifndef GUARD_lhe_io_hpp
#define GUARD_lhe_io_hpp

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>

// to use the existing structs
#include "showering.hpp"

// function to read momenta from an LHE particle line
Particle read_momenta(const std::string& inputline);

// read LHE file and extract particle momenta, event weights, and multiweights
void read_lhefile(const std::string& filename, 
                  std::vector<Event>& events,
                  std::vector<double>& weights, 
                  std::vector< std::map<std::string, double> >& multiweights);

// initialize LHE file for writing
std::ofstream init_lhe(const std::string& filename, double sigma, double stddev,
                    double ECM);

// write events to LHE file
void write_lhe(std::ofstream& outfile, const std::vector<Event>& events,
            double s_hat, bool debug = false);

// finalize LHE file
void finalize_lhe(std::ofstream& outfile);

#endif