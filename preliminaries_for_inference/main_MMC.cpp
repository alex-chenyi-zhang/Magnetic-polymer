#include "magpol.h"
#include <iostream>
#include <sstream>
#include <chrono>

using namespace std::chrono;
using namespace magpol;

int main(int argc, char* argv[]){
    // These all the input parameters that for now come from command line
    /*std::stringstream string1(argv[1]);
    std::stringstream string2(argv[2]);
    std::stringstream string3(argv[3]);
    std::stringstream string4(argv[4]);
    std::stringstream string5(argv[5]);
    std::stringstream string6(argv[6]);
    std::stringstream string7(argv[7]);

    int N_monomers;
    int stride_length;
    int number_of_strides;
    float J;
    float alpha;
    float inv_kbT;
    int n_temps;
    
    string1 >> N_monomers;
    string2 >> stride_length;
    string3 >> number_of_strides;
    string4 >> J;
    string5 >> alpha;
    string6 >> inv_kbT;
    string7 >> n_temps;


    std::cout << argv[0] << "\n";
    std::cout << "The number of monomers is: " << N_monomers << "\n";
    std::cout << "The length of the simulation is: " << number_of_strides*stride_length << "\n";
    std::cout << "The coupling values is: " << J << "\n";
    std::cout << "The value of alpha is: " << alpha << "\n";
    std::cout << "The value of the inverse temperature is: " << inv_kbT << "\n";
    std::cout << "The number of temperatures is: " << n_temps << "\n";
    




    multiple_markov_chains MMC_simulation(N_monomers, stride_length, number_of_strides, J, alpha, inv_kbT, n_temps);
    
    auto start = high_resolution_clock::now();
    
    MMC_simulation.run_MMC();
      
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << "\n";*/

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