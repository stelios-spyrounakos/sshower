#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>

#include <zlib.h>

#include "lhe_io.hpp"
// to use the existing structs
#include "showering.hpp"

using namespace std;

// read the momenta from the LHE file particle info line
Particle read_momenta(const string& inputline) {
    istringstream iss(inputline);
    Particle p;

    // split particle line into "tokens"
    /*
    vector<string> tokens;
    string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    */
    vector<string> tokens((istream_iterator<string>(iss)),
                        istream_iterator<string>());
    p.id = stoi(tokens[0]);
    p.status = stoi(tokens[1]);
    p.px = stod(tokens[6]);
    p.py = stod(tokens[7]);
    p.pz = stod(tokens[8]);
    p.E = stod(tokens[9]);
    p.m = stod(tokens[10]);
    // setting shower specific parameters for hard subprocess particles
    p.t_at_em = p.E * p.E;
    p.z_at_em = 1.0;
    p.stopped_evolving = false;
    
    return p;
}

// read LHE file and extract particle momenta, event weights, and multiweights
void read_lhefile(const string& filename, 
                  vector<Event>& events,
                  vector<double>& weights, 
                  vector< map<string, double> >& multiweights) {
    cout << "Opening LHE file " << filename << " for reading." << endl;
    // handle .gz compressed lhe files
    if (filename.substr(filename.size() - 3) == ".gz") {
        gzFile infile = gzopen(filename.c_str(), "rb");
        if (!infile) {
            cerr << "Failed to open the file: " << filename << endl;
            return;
        }
        char buffer[1024];
        bool reading_event = false;
        vector<Particle> ev_particles;
        double weight = 1.0;
        map<string, double> multiweight;
        while (gzgets(infile, buffer, sizeof(buffer))) {
            string line(buffer);

            if (line.find("<event>") != string::npos) {
                reading_event = true;
                ev_particles.clear();
                weight = 1.0;
                multiweight.clear();
            } 
            if (reading_event == true) {
                if (line.find("</event>") != string::npos) {
                    reading_event = false;
                    Event event {ev_particles};
                    events.push_back(event);
                    weights.push_back(weight);
                    multiweights.push_back(multiweight);
                }
                // split line into "tokens"
                istringstream iss(line);
                /*
                vector<string> tokens;
                string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }
                */
                vector<string> tokens((istream_iterator<string>(iss)),
                                    istream_iterator<string>());
                // extract the weight
                if (tokens.size() == 6) {
                    weight = stod(tokens[2]);
                }
                // extract the particle data
                if (tokens.size() == 13) {
                    ev_particles.push_back(read_momenta(line));
                }
                // extract multiweights
                if (tokens.size() == 4) {
                    string id = tokens[1];
                    size_t pos = id.find("id=");
                    id.erase(pos, 3);
                    id.erase(remove(id.begin(), id.end(), '>'), id.end());
                    id.erase(remove(id.begin(), id.end(), '\''), id.end());
                    multiweight[id] = stod(tokens[2]);
                }
            }
        }
        cout << "LHE file read." << endl;
        gzclose(infile);
    // handle uncompressed lhe files
    } else {
        ifstream infile(filename);
        if (!infile.is_open()) {
            cerr << "Failed to open the file: " << filename << endl;
            return;
        }
        string line;
        bool reading_event = false;
        vector<Particle> ev_particles;
        double weight = 1.0;
        map<string, double> multiweight;
        while (getline(infile, line)) {
            if (line.find("<event>") != string::npos) {
                reading_event = true;
                ev_particles.clear();
                weight = 1.0;
                multiweight.clear();
            } 
            if (reading_event == true) {
                if (line.find("</event>") != string::npos) {
                    reading_event = false;
                    Event event {ev_particles};
                    events.push_back(event);
                    weights.push_back(weight);
                    multiweights.push_back(multiweight);
                }
                // split line into "tokens"
                istringstream iss(line);
                /*
                vector<string> tokens;
                string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }
                */
                vector<string> tokens((istream_iterator<string>(iss)),
                                    istream_iterator<string>());
                // extract the weight
                if (tokens.size() == 6) {
                    weight = stod(tokens[2]);
                }
                // extract the particle data
                if (tokens.size() == 13) {
                    ev_particles.push_back(read_momenta(line));
                }
                // extract multiweights
                if (tokens.size() == 4) {
                    string id = tokens[1];
                    size_t pos = id.find("id=");
                    id.erase(pos, 3);
                    id.erase(remove(id.begin(), id.end(), '>'), id.end());
                    id.erase(remove(id.begin(), id.end(), '\''), id.end());
                    multiweight[id] = stod(tokens[2]);
                }
            }
        }
        cout << "LHE file read." << endl;
        infile.close();
    }
}

// initialize LHE file for writing
ofstream init_lhe(const string& filename, double sigma, double stddev,
                double ECM) {
    cout << "Opening LHE file " << filename << " for writing." << endl;
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Failed to open the file: " << filename << endl;
        return outfile;
    }
    // write header
    outfile << "<LesHouchesEvents version=\"1.0\">\n";
    outfile << "<!--\nFile generated with LHE C++ writer\n-->\n";
    outfile << "<init>\n";
    outfile << "\t11\t -11\t" << ECM / 2.0 << "\t" << ECM / 2.0 << "\t0\t0\t7\t7\t1\t1\n";
    outfile << "\t" << sigma << "\t" << stddev << "\t1.00000\t9999\n";
    outfile << "</init>\n";
    
    return outfile;
}

// write events to LHE file, given some particle information
// (information that was not and will not be used for our purposes will be set
// to the MadGraph lhe file value or a default that will be mentioned).
// The control flow structure checks are specific to e+e- -> qqbar
void write_lhe(ofstream& outfile, const vector<Event>& events,
            double s_hat, bool debug) {
    // loop through events
    for (size_t eventno = 0; eventno < events.size(); ++eventno) {
        Event event = events[eventno];

        // variables to store info for each event
        int ng = 0;
        vector<int> status;
        vector< vector<double> > momenta;
        vector<int> flavours;
        vector<int> colours;
        vector<int> anticolours;
        vector<double> helicities;
        vector< vector<int> > relations;

        // loop through the particles of the event
        for (size_t partno = 0; partno < event.particles.size(); ++partno) {
            Particle p = event.particles[partno];
            momenta.push_back({p.px, p.py, p.pz, p.E});
            status.push_back(p.status);
            // assign relations (same as MadGraph lhe)
            if (p.status == -1) {
                relations.push_back({0, 0});
            } else if (p.status == 1) {
                relations.push_back({1, 2});
            }
            // asssign flavours and helicities (helicity set to 1 for all)
            flavours.push_back(p.id);
            helicities.push_back(1);
            // asssign colours and anticolours
            if (abs(p.id) == 11) {
                colours.push_back(0);
                anticolours.push_back(0);
            } else if (abs(p.id) > 0 && abs(p.id) < 6) {
                if (p.id < 0) {
                    colours.push_back(0);
                    anticolours.push_back(501);
                } else {
                    colours.push_back(501);
                    anticolours.push_back(0);
                }
            } else if (p.id == 21) { // (gluons are singlets here)
                ng++;
                colours.push_back(500 + 2 * ng);
                anticolours.push_back(500 + 2 * ng);
            }
        }

        if (debug) {
            cout << "Writing event " << eventno + 1 << "\n";
            cout << "s_hat: " << s_hat << " (sqrt(s_hat): "
                << sqrt(s_hat) << ")\n";
        }

        // set form and precision
        outfile << scientific << setprecision(10);
        
        // write event header
        outfile << "<event>\n";
        outfile << momenta.size() << "\t9999\t1.000000\t" << sqrt(s_hat)
                << "\t0.0078125\t0.1187\n";

        // write particle information
        for (size_t i = 0; i < momenta.size(); ++i) {
            vector<double> momentum = momenta[i];
            vector<int> relation = relations[i];
            double mass = 0;  // we have set everything massless
            outfile << "\t" << flavours[i] << "\t" << status[i] << "\t"
                    << relation[0] << "\t" << relation[1] << "\t" << colours[i]
                    << "\t" << anticolours[i] << "\t" << momentum[0] << "\t"
                    << momentum[1] << "\t" << momentum[2] << "\t" << momentum[3]
                    << "\t" << mass << "\t0.0000e+00\t" << helicities[i] << "\n";
            if (debug) {
                cout << "\t" << flavours[i] << "\t" << status[i] << "\t"
                    << relation[0] << "\t" << relation[1] << "\t" << colours[i]
                    << "\t"  << anticolours[i] << "\t" << momentum[0] << "\t"
                    << momentum[1] << "\t" << momentum[2] << "\t" << momentum[3]
                    << "\t" << mass << "\t0.0000e+00\t" << helicities[i]
                    << "\n";
            }
        }
        // write closing event tag
        outfile << "</event>\n";
    }
}

// finalize LHE file
void finalize_lhe(ofstream& outfile) {
    outfile << "</LesHouchesEvents>\n";
    outfile.close();
    cout << "LHE file closed." << endl;
}