#include <iostream>
#include <cmath>
#include <iomanip>

#include "progress_bar.hpp"

using namespace std;

// function to display the progress bar (input is a number from 0 to 1)
void show_progress_bar(double progress) {
    static int last_percentage = -1;
    int percentage = round(progress * 100);
    
    if (percentage != last_percentage) {
        int bar_width = 70;
        int prog_position = round(bar_width * progress);
        
        cout << "[";
        // display the filled part
        for (int i = 0; i < bar_width; ++i) {
            if (i < prog_position) { cout << "#"; }
            else if (i == prog_position) { cout << ".."; }
            else { cout << " "; }
        }
        // \r to overwrite the same line
        cout << "] " << setw(3) << percentage << " %\r";
        // flush to write immediately without waiting for the buffer to fill
        cout.flush();
        last_percentage = percentage;
    }
    if (progress >= 1.0) {
        cout << endl;
        last_percentage = -1;
    }
}
