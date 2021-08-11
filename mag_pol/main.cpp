#include "magpol.h"
#include <iostream>
#include <sstream>
#include <chrono>

using namespace std::chrono;
using namespace magpol;

int main(int argc, char* argv[]){
//int main(){
    std::stringstream string1(argv[1]);
    std::stringstream string2(argv[2]);
    std::stringstream string3(argv[3]);
    int N_monomers;
    int simul_length;
    float J;
    string1 >> N_monomers;
    string2 >> simul_length;
    string3 >> J;
    std::cout << argv[0] << "\n";
    std::cout << "The number of monomers is: " << N_monomers << "\n";
    std::cout << "The length of the simulation is: " << simul_length << "\n";
    std::cout << "The coupling values is: " << J << "\n";
    
    auto start = high_resolution_clock::now();
    saw_MC simulation(N_monomers, simul_length, J); 
    simulation.run();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << "\n";

    //simulation.write_results_on_file();
    return 0;
}