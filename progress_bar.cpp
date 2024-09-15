#include <iostream>
#include <chrono>
#include <thread>

using namespace std;

// function to display the progress bar (input is a number from 0 to 1)
void show_progress_bar(float progress) {
    int bar_width = 70;
    int prog_position = bar_width * progress;
    int percentage = progress * 100;
    
    cout << "[";
    // display the filled part
    for (int i = 0; i < bar_width; ++i) {
        if (i < prog_position) { cout << "="; }
        else if (i == prog_position) { cout << ">"; }
        else { cout << " "; }
    }
    // \r to overwrite the same line
    cout << "] " << percentage << " %\r";
    // flush to write immediately without waiting for the buffer to fill
    cout.flush();
}
