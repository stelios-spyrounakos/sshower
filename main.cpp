#include <cmath>
#include <iostream>
#include <vector>
//#include <algorithm>
//#include <iterator>
//#include <ctime>

#include "lhe_io.hpp"
#include "showering.hpp"
#include "kinematic_reco.hpp"
#include "progress_bar.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cerr << "An input filename and an output filename need to be "
            << "provided as arguments" << endl;
        return 1;
    }
    else if (argc > 3) {
        cerr << "Only input filename and output filename are needed "
            << "The rest oof the arguments are not used." << endl;
    }
    string infilename = argv[1];
    string outfilename = argv[2];
    cout << "Showering " << infilename << endl;

    // read input lhe file
    vector<Event> events;
    vector<double> weights;
    vector< map<string, double> > multiweights;
    read_lhefile(infilename, events, weights, multiweights);
    int n_events = 100; // events.size();

    // shower events
    double Q_cutoff = 0.935;    // (to match with Herwig 7)
    vector<Event> showered_events;
    bool check_mom_conservation = true;
    showered_events.reserve(n_events);
    for (int i = 0; i < n_events; ++i) {
        Event showered_event;

        ///// FOR DEBUGGING /////
        cout << "SHOWERING EVENT " << i << endl;
        /////////////////////////

        showered_event.particles = shower_event(events[i], Q_cutoff);

        ///// FOR DEBUGGING /////
        cout << "SHOWERED EVENT " << i << endl;
        /////////////////////////

        showered_events.push_back(showered_event);

        // if needed, checks if the 3-momentum of the event is conserved
        vector<double> total_mom(3);
        total_mom = check_mom_cons(showered_event);
        if (check_mom_conservation == true) {
            cout << "Total momentum of the event after showering: ("
                << total_mom[0] << ", " << total_mom[1] << ", " << total_mom[3]
                << ")." << endl;
        }

        show_progress_bar(i / (n_events - 1));

        ///// FOR DEBUGGING /////
        cout << endl << endl << endl << endl << endl;
        /////////////////////////
    }

    // write output lhe file
    double sigma = 1.2;
    double error = 0.2;
    double E_CM = 206;
    double s_hat = E_CM * E_CM;
    ofstream outfile = init_lhe(outfilename, sigma, error, E_CM);
    write_lhe(outfile, showered_events, s_hat);
    finalize_lhe(outfile);

    return 0;
}