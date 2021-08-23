#include "magpol.h"
#include <iostream>
#include <sstream>
#include <chrono>

using namespace std::chrono;
using namespace magpol;

int main(int argc, char* argv[]){
    // These all the input parameters that for now come from command line
    std::stringstream string1(argv[1]);
    std::stringstream string2(argv[2]);
    std::stringstream string3(argv[3]);
    std::stringstream string4(argv[4]);
    std::stringstream string5(argv[5]);
    std::stringstream string6(argv[6]);
    std::stringstream string7(argv[7]);

    int n_temps;// = 20;                   // number of temperatures
    std::random_device rand_dev;                         // generates a random seed, to initialize a random engine
    std::mt19937 mer_twist;
    int N_monomers;
    int stride_length;
    int number_of_strides;
    float J;
    float alpha;
    float inv_kbT;
    float max_inv_temp = 0.3;
    
    string1 >> N_monomers;
    string2 >> stride_length;
    string3 >> number_of_strides;
    string4 >> J;
    string5 >> alpha;
    string6 >> inv_kbT;
    string7 >> n_temps;

    std::uniform_real_distribution<double> uni_01(0.0, 1.0);  // uniform in [0,1]. Used to acc/rej MC move
    std::uniform_int_distribution<int> uni_temps(0, n_temps-2);


    int * temporary_spins1 = new int[N_monomers];
    int * temporary_spins2 = new int[N_monomers];
    int ** temporary_coords1 = new int*[N_monomers];
    int ** temporary_coords2 = new int*[N_monomers];
    for (int i_mono = 0; i_mono < N_monomers; i_mono++){
        temporary_coords1[i_mono] = new int[3];
        temporary_coords2[i_mono] = new int[3];
    }

    /*for (int i_mono = 0; i_mono < N_monomers; i_mono++){
        for (int j = 0; j < 3; j++){
            std::cout << temporary_coords1[i_mono][j] << " ";
        }
        std::cout << "\n";
    }*/

    std::cout << argv[0] << "\n";
    std::cout << "The number of monomers is: " << N_monomers << "\n";
    std::cout << "The length of the simulation is: " << number_of_strides*stride_length << "\n";
    std::cout << "The coupling values is: " << J << "\n";
    std::cout << "The value of alpha is: " << alpha << "\n";
    std::cout << "The value of the inverse temperature is: " << inv_kbT << "\n";
    std::cout << "The number of temperatures is: " << n_temps << "\n";
    
    
    auto start = high_resolution_clock::now();
    saw_MC** simulations = new saw_MC*[n_temps];
    // do many parallel simulations with inverse temperatures uniformly spaced between 0 and the value of interest
    for (int i = 0; i < n_temps; i++){
        simulations[i] = new saw_MC(N_monomers, stride_length, number_of_strides, J, alpha, inv_kbT-(inv_kbT-max_inv_temp)*i/(n_temps-1)); 
    }
    int accepted_swaps = 0;

    for (int i_strides = 0; i_strides < number_of_strides; i_strides++){
        for (int i_temps = 0; i_temps < n_temps; i_temps++){
            simulations[i_temps]->run();
        }

        /* now do the actual MMC monte carlo part where you pick a number n between 0 and n_temps-2
           and try to swap the states of the n-th and the (n+1)-th markov chain*/
        int swap = uni_temps(mer_twist);
        float alpha_swap;
        //std::cout << "chain to be swapped:  " << chain_to_swap << "\n";
        float delta_en = 0;
        delta_en += simulations[swap]->get_beta() * (simulations[swap+1]->get_ene() - simulations[swap]->get_ene());
        delta_en += simulations[swap+1]->get_beta() * (simulations[swap]->get_ene() - simulations[swap+1]->get_ene());
        if (delta_en <= 0){ alpha_swap = 1; }
        else{ alpha_swap = exp(-delta_en); }

        if (alpha_swap >= uni_01(mer_twist)){
            simulations[swap]->copy_spins(temporary_spins1);
            simulations[swap+1]->copy_spins(temporary_spins2);
            simulations[swap]->copy_coords(temporary_coords1);
            simulations[swap+1]->copy_coords(temporary_coords2);

            /*for (int i_mono = 0; i_mono < N_monomers; i_mono++){
                for (int j = 0; j < 3; j++){
                    std::cout << temporary_coords1[i_mono][j] << " ";
                }
                std::cout << "\n";
            }*/

            simulations[swap]->set_spins(temporary_spins2);
            simulations[swap+1]->set_spins(temporary_spins1);
            simulations[swap]->set_coords(temporary_coords2);
            simulations[swap+1]->set_coords(temporary_coords1);
            
            accepted_swaps++;
        }
    }

    std::cout << "The number of accepted swaps was " << accepted_swaps << " out of " << number_of_strides << "\n";

    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << "\n";

    for (int i_temps = 0; i_temps < n_temps; i_temps++){
        simulations[i_temps]->write_results_on_file();
    }






    /*Test** arr = new Test*[5];
    for (int i = 0; i < 5; i++){
        arr[i] = new Test(i, i+1);
    }
    for (int i = 0; i < 5; i++){
        arr[i]->print();
    }*/
    for (int i_mono = 0; i_mono < N_monomers; i_mono++){
        delete [] temporary_coords1[i_mono];
        delete [] temporary_coords2[i_mono];
    }
    delete [] temporary_spins1;
    delete [] temporary_spins2;
    delete [] temporary_coords1;
    delete [] temporary_coords2;

    return 0;
}