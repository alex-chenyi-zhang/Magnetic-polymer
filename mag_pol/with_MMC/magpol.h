#ifndef MAGPOL_H
#define MAGPOL_H

#include <random>
namespace magpol{
class linear_hash{
    public:
        int M = 5000000;
        int a[3] = {17,290,4914};
        int * occupancy;
        int * monomer_index;
        int ** key_values;
        linear_hash();
        ~linear_hash();
};

class saw_MC{
    int stride, n_swaps;                           // in the Multiple Markov Chain stride is the number of MC moves before attempting a swap
                                                   // n_steps is set equal to n_swaps*stride 
    int done_strides = 0;                          // incremented every time a stride is done. Needed to properly save the observables we need
    int n_mono, n_steps;
    int ** coord;                                  // 2D array that contains the coordinates of the monomers
    int ** trial_coord;
    int * Ree2;                                    // Here we store the squared end-to-end distance
    float *Rg2;
    int * spins;                                   // Here I store the spin congiguration of the Ising/Potts subsystem
    int * trial_spins;
    float * h_fields;                              // Values of the local fields in the Ising/Potts model
    float max_field_value = 2;
    float spin_coupling;
    float energy = 0;
    float alpha_h;                           // alpha is between 0 and 1. it is the importance of the polymer part of the hamiltonian vs the 
                                                   // the part that ONLY depends on the spins i.e. the parts with the external fields
    float beta_temp;
    float * energies;
    float * magnetization;
    int p_moves[47][2][3] = {};                // The 47 orthogonal transformations of the octahedral group (except identity)

    int ** neighbours;                             /* This is the list of the neighbours of each monomer. Needed to evaluate the hamiltonian.          
                                                      To compute energy is order n_mono, so is the same order of the other operations needed
                                                      to propose a move, so I guess it's ok. To increase the efficiency one should compute only
                                                      the delta energy, but it's not so straighforward as the pivot is a global move. In this 
                                                      probably a trial_neighbour matrix will be needed*/
    int ** trial_neighbours;

    std::random_device rd;                         // generates a random seed, to initialize a random engine
    std::mt19937 mt;
    std::uniform_real_distribution<double> uni_R;  // uniform in [0,1]. Used to acc/rej MC move
    std::uniform_int_distribution<int> uni_I_poly; // uniform in [1,n_mono-1]. To pick a random monomer. We don't do pivots around the 0-th monomer 
    std::uniform_int_distribution<int> uni_I_poly2;// uniform in [1,n_mono-3]. Needed for the one-bead-flip
    std::uniform_int_distribution<int> uni_I_poly3;// uniform in [1,n_mono-4]. Needed for cranckshaft-type moves
    std::uniform_int_distribution<int> uni_G;      // uniform in [0,47-1]. To pick a random element of the symmetry group
    std::uniform_int_distribution<int> uni_spins;  // random number picked between {-1,0,+1}
    std::uniform_int_distribution<int> uni_spins2; // random number picked between {0,1}
    std::uniform_int_distribution<int> local_move_rand;  // allows you to choose one of the 4 implemented local moves
    int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int tr_signs[8][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},       
    {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};      // by combining perm's and sign comb's you obtain the 48 pivot moves
    linear_hash hash_saw;                          // We declare here the hash table that will be used for self-avoidance checks
    int * hashed_where;                            // at each attempted pivot you store here the hashed coordinates of the monomers
                                                   // :this is useful for quick cleanup of the hashtable.
    int * whos_hashed;                             // constains the sequence of monomers used for a self-av check
    int n_hashes;                                  // Tells you in each self_av. check how many monomers you inserted in the hash table
    public:
        saw_MC(int,int,int,float,float,float);                           //constructor
        ~saw_MC();                                 //destructor
        void try_pivot(int,int);
        int hash_function(int*,int,int*);
        bool check_saw(int); 
        void compute_neighbour_list(int**, int**);          // if argument is 1 uses coord, if it is 0 uses trial_coord  
        float compute_new_energy(int**);
        void spins_MC(void);
        void remove_from_hash_table(int);
        bool add_to_hash_table(int, int*);
        bool check_site_occupancy(int, int*);
        void single_bead_flip(void);
        void crankshaft_180(void);
        void crankshaft_90_270(int);
        void run(void);      
        void write_results_on_file(void);   
        float gyr_rad_square(void);         
};



}
#endif // MAGPOL_H