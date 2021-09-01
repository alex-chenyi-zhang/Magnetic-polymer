#include "magpol.h"
#include <iostream>
#include <sstream>
#include <chrono>

using namespace std::chrono;
using namespace magpol;

int main(int argc, char* argv[]){
    std::stringstream string1(argv[1]);
    std::string input_filename;
    string1 >> input_filename;

    multiple_markov_chains MMC_simulation(input_filename);
    
    auto start = high_resolution_clock::now();
    
    MMC_simulation.run_MMC();
      
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << "\n";


    return 0;
}