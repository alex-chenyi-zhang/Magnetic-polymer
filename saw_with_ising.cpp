#include <iostream>
#include <random>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

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



linear_hash::linear_hash(){
    occupancy = new int[M];
    monomer_index = new int[M];
    key_values = new int*[M];
    for (int i = 0; i < M; i++){
        key_values[i] = new int[3];
        occupancy[i] = 0;
        monomer_index[i] = -1;
        for (int j = 0; j < 3; j++){
            key_values[i][j] = 0;
        }
    }
}



linear_hash::~linear_hash(){
    for (int i = 0; i < M; i++){
        delete [] key_values[i];
    }
    delete [] key_values;
    delete [] occupancy;
    delete [] monomer_index;
}



class saw_MC{
    int n_mono, n_steps;
    int ** coord;                                  // 2D array that contains the coordinates of the monomers
    int ** trial_coord;
    int * Ree2;                                    // Here we store the squared end-to-end distance
    int * spins;                                   // Here I store the spin congiguration of the Ising/Potts subsystem
    int * trial_spins;
    float * h_fields;                              // Values of the local fields in the Ising/Potts model
    float max_field_value = 2;
    float spin_coupling;
    float energy = 0;
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
        saw_MC(int,int,float);                           //constructor
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
};



/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

/* The constructor is called in the beginning and initializes the
   attributes of the class. It takes for now as parameters the 
   number of MC steps the simulation lasts and the length of the polymer
*/
saw_MC::saw_MC(int x, int y, float j_spins) : uni_R(0.0, 1.0), uni_I_poly(1,x-1), uni_G(0,47-1), uni_spins(-1,1), 
                                uni_I_poly2(0,x-3), uni_I_poly3(0,x-4), local_move_rand(0,3), uni_spins2(0,1), mt(1){   //CONSTRUCTOR // mt(rd())--> This for a random seed
    n_mono = x;
    n_steps = y; 
    spin_coupling = j_spins;
    coord = new int*[n_mono];
    trial_coord = new int*[n_mono];
    Ree2 = new int[n_steps];
    energies = new float[n_steps];
    magnetization = new float[n_steps];
    spins = new int[n_mono];
    trial_spins = new int[n_mono];
    h_fields = new float[n_mono];
    neighbours = new int*[n_mono];
    trial_neighbours = new int*[n_mono];
    hashed_where = new int[n_mono];
    whos_hashed = new int[n_mono];
    for (int i = 0; i < n_mono; i++){
        coord[i] = new int[3];
        trial_coord[i] = new int[3];
        neighbours[i] = new int[7]; // 1st entry: number of neighbours--> the others are the neighbours if there are any
        trial_neighbours[i] = new int[7];
    }
    for (int i = 0; i < n_mono; i++){
        coord[i][0] = i;
        coord[i][1] = 0;
        coord[i][2] = 0;
        trial_coord[i][0] = i;
        trial_coord[i][1] = 0;
        trial_coord[i][2] = 0;
        spins[i] = uni_spins(mt);//*2-1;
        trial_spins[i] = spins[i];
        h_fields[i] = 0;/*(uni_R(mt)-1)*max_field_value*2;  LET'S USE O-FIELDS FOR NOW TO SEE IF WE CAN QUALITATIVELY GET THE RESULTS FROM COLI' ET AL*/
        neighbours[i][0] = 0;     // starting from a straight rod no one has any neighbour
        trial_neighbours[i][0] = 0;
        for (int j = 0; j < 5; j++){
            neighbours[i][j+1] = -1; // -1 means it's not associated to any of the mono's that are indexed from 0 to n_mono-1
            trial_neighbours[i][j+1] = -1;
        }
    } // I initialize the structure as a fully extended polymer in the x-direction
      // with the first monomer anchored in the origin.

      // MAGARI A UN CERTO PUNTO PROVA A VEDERE LA DIMERIZATION COME INIZIALIZZAZIONE (VEDI SOKAL E MADRAS)

    std::cout << "You have initialized your polymer configuration!!!\n";

    // Now initialize the array containing the information for the different transformations
    int count = 0;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 8; j++){
            if(i==0 && j==0) {continue;} //we don't care about identity
            for (int k = 0; k < 3; k++){
                //p_moves[count][k][perms[i][k]] = tr_signs[j][k];
                p_moves[count][0][k] = perms[i][k];
                p_moves[count][1][k] = tr_signs[j][k];
            }
            count++;
        }
    }
}


saw_MC::~saw_MC(){    // DESTRUCOR
    for (int i = 0; i < n_mono; i++){
        delete [] coord[i];
        delete [] trial_coord[i];
        delete [] neighbours[i];
        delete [] trial_neighbours[i];
    }
    delete [] coord;
    delete [] trial_coord;
    delete [] Ree2;
    delete [] energies;
    delete [] magnetization;
    delete [] hashed_where;
    delete [] whos_hashed;
    delete [] spins;
    delete [] trial_spins;
    delete [] h_fields;
    delete [] neighbours;
    delete [] trial_neighbours;
}



/*********************************************************************************************/
/* This function makes the g-th transformation on the k-th pivot point
   and stores the transformed configuration in trial_coord. One then 
   needs to check the self-avoidance of this new walk. If it OK it 
   overwrites the old walk stored in coord.
*/
void saw_MC::try_pivot(int k, int g){
    int pivot = k;
    int symm  = g; 

    for (int i = 0; i < pivot+1; i++){   // The monomers before the pivot point do not move
        for (int j = 0; j < 3; j++){
            trial_coord[i][j] = coord[i][j];
        }
    }

    for (int i = pivot + 1; i < n_mono; i++){
        for (int j = 0; j < 3; j++){
            trial_coord[i][j] = p_moves[symm][1][j]*(coord[i][p_moves[symm][0][j]] - coord[pivot][p_moves[symm][0][j]]) + coord[pivot][j];
        }
    }
}



/*********************************************************************************************/
int saw_MC::hash_function(int a[], int M, int coord_to_be_hashed[]){
    int sum = 0;
    for (int i = 0; i < 3; i++){
        sum += a[i]*coord_to_be_hashed[i];
    }
    return abs(sum)%M;
}




/*********************************************************************************************/
bool saw_MC::check_site_occupancy(int i_mono, int coord_to_hash[]){
    bool already_there = false;
    int h_coord = hash_function(hash_saw.a, hash_saw.M, coord_to_hash);
    for (int i = h_coord; i < hash_saw.M+h_coord; i++){
        if(hash_saw.occupancy[i%hash_saw.M] == 0){ 
            break;
        }
        else{
            already_there = true;
            for (int j = 0; j< 3; j++){
                already_there = already_there && (hash_saw.key_values[i%hash_saw.M][j] == coord_to_hash[j]);
            }
            if (already_there){
                break;
            }
        }
    }
    return already_there;

}
/*********************************************************************************************/
// This adds a lattice site to the hash table and tells you whether the site was already inserted or not
bool saw_MC::add_to_hash_table(int i_mono, int coord_to_hash[]){
    bool already_there = false;

    int h_coord = hash_function(hash_saw.a, hash_saw.M, coord_to_hash);
    for (int i = h_coord; i < hash_saw.M+h_coord; i++){
        if(hash_saw.occupancy[i%hash_saw.M] == 0){ 
            hash_saw.occupancy[i%hash_saw.M] = 1;
            hash_saw.monomer_index[i%hash_saw.M] = i_mono;
            for (int j = 0; j< 3; j++){
                    hash_saw.key_values[i%hash_saw.M][j] = coord_to_hash[j];
            }
            hashed_where[i_mono] = i%hash_saw.M;
            break;
        }
        else{
            already_there = true;
            for (int j = 0; j< 3; j++){
                already_there = already_there && (hash_saw.key_values[i%hash_saw.M][j] == coord_to_hash[j]);
            }
            if (already_there){
                break;
            }
        }
    }
    return already_there;
}



/*********************************************************************************************/
void saw_MC::remove_from_hash_table(int hash_address){
    if(hash_saw.occupancy[hash_address] != 0){
        hash_saw.occupancy[hash_address] = 0;   // delete the element you actually want to delete;
        for (int i = hash_address + 1; i < hash_saw.M + hash_address + 1; i++){
            if(hash_saw.occupancy[i] != 0){     // If the next element in the hash table is occupied, delete the element and insert it from scratch;
                hash_saw.occupancy[i] = 0;
                add_to_hash_table(hash_saw.monomer_index[i], hash_saw.key_values[i]);
            }
            else{
                break;
            }
        }
    }
    else{
        std::cout << "ERROR!! The element you're trying to remove does not exist!! \n";
    }
}


/*********************************************************************************************/
/* The input is the pivot point used to generate the trial walk. 
   it is useful to know, as it is more convenient to check self
   avoidance from the pivot point outwards.
*/
bool saw_MC::check_saw(int k){
    int L = std::max(k,n_mono-1-k) + 1;
    bool is_saw = true;
    int t = 0;
    n_hashes = 0;
    bool there_already;
    while(t<L && is_saw){
        if (k+t < n_mono){
            there_already = add_to_hash_table(k+t, trial_coord[k+t]);
            if (there_already){ is_saw = false; }
            else{
                whos_hashed[n_hashes] = k+t;
                n_hashes++;
            }
        }
        if (k-t >= 0 && t != 0){
            there_already = add_to_hash_table(k-t, trial_coord[k-t]);
            if (there_already){ is_saw = false; }
            else{
                whos_hashed[n_hashes] = k-t;
                n_hashes++;
            }
        }

        t++;
    }
    return is_saw;
}



/*********************************************************************************************/
void saw_MC::single_bead_flip(){
    int mono = uni_I_poly2(mt);
    //std::cout << mono << "\n";
    int dist = 0;
    for (int j = 0; j < 3; j++){
        dist += (trial_coord[mono+2][j] - trial_coord[mono][j])*(trial_coord[mono+2][j] - trial_coord[mono][j]);
    }
    if (dist == 2){  // i.e. if the bead in between form a bend of the polymer
        int address_to_empty = hashed_where[mono+1];
        int new_coord[3];
        for (int j = 0; j < 3; j++){
            new_coord[j] = trial_coord[mono][j] + trial_coord[mono+2][j] - trial_coord[mono+1][j];
        }

        bool there_already;
        there_already = check_site_occupancy(mono+1, new_coord);
        if (!there_already){
            for (int j = 0; j < 3; j++){
                trial_coord[mono+1][j] = new_coord[j];
            }
            remove_from_hash_table(address_to_empty);
            add_to_hash_table(mono+1, new_coord);
        }
    }
    
}



/*********************************************************************************************/
void saw_MC::crankshaft_90_270(int verso){
    // verso +1 means clockwise (90 degrees), -1 is 270
    int mono = uni_I_poly3(mt);
    int dist = 0;
    for (int j = 0; j < 3; j++){
        dist += (trial_coord[mono+3][j] - trial_coord[mono][j])*(trial_coord[mono+3][j] - trial_coord[mono][j]);
    }
    if (dist == 1){  // i.e. the chosen 3 segment long subchain forms a U-shape
        int new_coord1[3];
        int new_coord2[3];
        int axis;
        int orient;
        for (int j = 0; j < 3; j++){
            if ((trial_coord[mono+3][j] - trial_coord[mono][j]) != 0){
                axis = j;
                orient = trial_coord[mono+3][j] - trial_coord[mono][j];
                break;
            }
        }

        switch (axis) {
            case 0:
                new_coord1[0] = trial_coord[mono+1][0];
                new_coord1[1] = trial_coord[mono][1] - (trial_coord[mono+1][2]-trial_coord[mono][2])*orient*verso;
                new_coord1[2] = trial_coord[mono][2] + (trial_coord[mono+1][1]-trial_coord[mono][1])*orient*verso;

                new_coord2[0] = trial_coord[mono+2][0];
                new_coord2[1] = trial_coord[mono+3][1] - (trial_coord[mono+2][2]-trial_coord[mono+3][2])*orient*verso;
                new_coord2[2] = trial_coord[mono+3][2] + (trial_coord[mono+2][1]-trial_coord[mono+3][1])*orient*verso;
                break;
            case 1:
                new_coord1[0] = trial_coord[mono][0] + (trial_coord[mono+1][2]-trial_coord[mono][2])*orient*verso;
                new_coord1[1] = trial_coord[mono+1][1];
                new_coord1[2] = trial_coord[mono][2] - (trial_coord[mono+1][0]-trial_coord[mono][0])*orient*verso;

                new_coord2[0] = trial_coord[mono+3][0] + (trial_coord[mono+2][2]-trial_coord[mono+3][2])*orient*verso;
                new_coord2[1] = trial_coord[mono+2][1];
                new_coord2[2] = trial_coord[mono+3][2] - (trial_coord[mono+2][0]-trial_coord[mono+3][0])*orient*verso;
                break;
            case 2:
                new_coord1[0] = trial_coord[mono][0] - (trial_coord[mono+1][1]-trial_coord[mono][1])*orient*verso;
                new_coord1[1] = trial_coord[mono][1] + (trial_coord[mono+1][0]-trial_coord[mono][0])*orient*verso;
                new_coord1[2] = trial_coord[mono+1][2];

                new_coord2[0] = trial_coord[mono+3][0] - (trial_coord[mono+2][1]-trial_coord[mono+3][1])*orient*verso;
                new_coord2[1] = trial_coord[mono+3][1] + (trial_coord[mono+2][0]-trial_coord[mono+3][0])*orient*verso;
                new_coord2[2] = trial_coord[mono+2][2];
                break;
            default:
                std::cout << "!!! Invalid rotation axis !!! \n";
        }

        bool there_already1;
        bool there_already2;
        there_already1 = check_site_occupancy(mono+1, new_coord1);
        there_already2 = check_site_occupancy(mono+2, new_coord2);
        if (!there_already1 && !there_already2){
            //std::cout << "WOOOOOOOW\n\n";
            for (int j = 0; j < 3; j++){
                trial_coord[mono+1][j] = new_coord1[j];
                trial_coord[mono+2][j] = new_coord2[j];
            }
            int address_to_empty1 = hashed_where[mono+1];
            remove_from_hash_table(address_to_empty1);
            int address_to_empty2 = hashed_where[mono+2];
            remove_from_hash_table(address_to_empty2);

            add_to_hash_table(mono+1, new_coord1);
            add_to_hash_table(mono+2, new_coord2);
        }
    }
}

/*********************************************************************************************/
void saw_MC::crankshaft_180(){
    int mono = uni_I_poly3(mt);
    int dist = 0;
    for (int j = 0; j < 3; j++){
        dist += (trial_coord[mono+3][j] - trial_coord[mono][j])*(trial_coord[mono+3][j] - trial_coord[mono][j]);
    }
    if (dist == 1){  // i.e. the chosen 3 segment long subchain forms a U-shape
        int new_coord1[3];
        int new_coord2[3];
        for (int j = 0; j < 3; j++){
            new_coord1[j] = 2 * trial_coord[mono][j] - trial_coord[mono+1][j];
            new_coord2[j] = 2 * trial_coord[mono+3][j] - trial_coord[mono+2][j];
        }

        bool there_already1;
        bool there_already2;
        there_already1 = check_site_occupancy(mono+1, new_coord1);
        there_already2 = check_site_occupancy(mono+2, new_coord2);
        if (!there_already1 && !there_already2){
            //std::cout << "WOOOOOOOW\n\n";
            for (int j = 0; j < 3; j++){
                trial_coord[mono+1][j] = new_coord1[j];
                trial_coord[mono+2][j] = new_coord2[j];
            }
            int address_to_empty1 = hashed_where[mono+1];
            remove_from_hash_table(address_to_empty1);
            int address_to_empty2 = hashed_where[mono+2];
            remove_from_hash_table(address_to_empty2);

            add_to_hash_table(mono+1, new_coord1);
            add_to_hash_table(mono+2, new_coord2);
        }
    }
}
/*********************************************************************************************/
void saw_MC::compute_neighbour_list(int **coo, int **near){ // this method is run only if the pivot was successful, otherwise it doesn't make sense
    for (int i_mono = 0; i_mono < n_mono; i_mono++){
        near[i_mono][0] = 0;  // number of neighbours of the i_mono-th monomer
        for (int i = 1; i < 7; i++) { near[i_mono][i] = -1; } // i starts from 1 here
        
        int neigh_coordinates[3];
        for (int i = 0; i < 3; i++){ neigh_coordinates[i] = coo[i_mono][i];}
        
        
        for (int j = 0; j < 3; j++){
            for (int k = -1; k <= 1; k=k+2){
                neigh_coordinates[j] += k;

                int hash_neigh = hash_function(hash_saw.a, hash_saw.M, neigh_coordinates);
                for (int l = hash_neigh; l < hash_saw.M + hash_neigh; l++){
                    if(hash_saw.occupancy[l%hash_saw.M] == 0){break;}
                    else{
                        bool is_neighbour = true;
                        for (int m = 0; m< 3; m++){
                            is_neighbour = is_neighbour && (hash_saw.key_values[l%hash_saw.M][m] == neigh_coordinates[m]);
                        }
                        if (is_neighbour){
                            near[i_mono][0] += 1;
                            int n_neigh = near[i_mono][0];
                            near[i_mono][n_neigh] = hash_saw.monomer_index[l%hash_saw.M];
                            break;
                        }
                    }

                }

                neigh_coordinates[j] -= k;
            }
        }
    }
}

float saw_MC::compute_new_energy(int **near){
    float ENE = 0;
    for (int i_mono = 0; i_mono < n_mono; i_mono++){
        ENE = ENE - h_fields[i_mono] * spins[i_mono];
        for (int j = 0; j < near[i_mono][0]; j++){
            if (( spins[i_mono] == spins[near[i_mono][j+1]] ) && (spins[i_mono] != 0)) { 
                ENE = ENE - spin_coupling/2;
            }
        }
    }
    return ENE;
}


void saw_MC::spins_MC(){
    int n_flips = n_mono;

    // Now you can do as many spin flips as you wish with an annealed polymer configuration
    float acc;
    float delta_ene;
    float local_ene;
    float trial_local_ene;
    int trial_spin_value;
    int flip_candidate;
    for(int i_flip = 0; i_flip < n_flips; i_flip++){
        flip_candidate = uni_I_poly(mt); // pick a random spin to try to flip
        
        local_ene = 0;
        local_ene = local_ene - h_fields[flip_candidate] * spins[flip_candidate];
        for (int j = 0; j < neighbours[flip_candidate][0]; j++){
            if (( spins[flip_candidate] == spins[neighbours[flip_candidate][j+1]] ) && (spins[flip_candidate] != 0)) { 
                local_ene = local_ene - spin_coupling;
            }
        }
    
        if(spins[flip_candidate] == -1){ trial_spin_value = uni_spins2(mt); } //pick between 0 and 1
        else if(spins[flip_candidate] == 1) { trial_spin_value = -uni_spins2(mt); } //-------0 and -1
        else { trial_spin_value = uni_spins2(mt)*2 - 1; }                     //------------ 1 and -1

        trial_local_ene = 0;
        trial_local_ene = trial_local_ene - h_fields[flip_candidate] * trial_spin_value;
        for (int j = 0; j < neighbours[flip_candidate][0]; j++){
            if (( trial_spin_value == spins[neighbours[flip_candidate][j+1]] ) && (trial_spin_value != 0)) { 
                trial_local_ene = trial_local_ene - spin_coupling;
            }
        }
        delta_ene = trial_local_ene - local_ene;
        //std::cout << delta_ene << "\n";

        if (delta_ene <= 0){ acc = 1; }
        else {acc = exp(-delta_ene); }

        if (acc > uni_R(mt)){
            energy = energy + delta_ene;
            spins[flip_candidate] = trial_spin_value;
        }
    }
}


void saw_MC::run(){
    std::cout << check_saw(1) << "\n";
    compute_neighbour_list(coord, neighbours);
    for (int i = 0; i < n_hashes; i++){
        hash_saw.occupancy[hashed_where[whos_hashed[i]]] = 0;
    }
    energy = compute_new_energy(neighbours);
    //std::cout << energy << "\n";

    int n_acc = 0;
    int n_pivots = 0;
    int pivot_point;
    int transformation;
    bool is_still_saw;
    float trial_energy;
    float delta_energy;
    float acceptance;
    float magnet;
    for (int i = 0; i < n_steps; i++){
        Ree2[i] = coord[n_mono-1][0]*coord[n_mono-1][0]+coord[n_mono-1][1]*coord[n_mono-1][1]+coord[n_mono-1][2]*coord[n_mono-1][2];
        energies[i] = energy;
        magnet = 0;
        for (int i_mono = 0; i_mono < n_mono; i_mono++){
            magnet += spins[i_mono];
        }
        magnetization[i] = magnet/n_mono;
        //std::cout << energy -compute_new_energy() << "\n";
        pivot_point    = uni_I_poly(mt);
        transformation = uni_G(mt);
        try_pivot(pivot_point,transformation);
        is_still_saw = check_saw(pivot_point);
        if (i%10000 == 0){
            std::cout << "i_step = "<< i << ' '<<pivot_point << ' ' << transformation <<' '<< is_still_saw << "\n";
        }

        /* In this next section if the pivot move yielded an invalid saw, then in order to do 
           the local moves efficiently we insert the pre-pivot config. in the hash table*/
        if (!is_still_saw){ // It the new config is not saw we try the local moves on the previous config.
            for (int i = 0; i < n_hashes; i++){
                hash_saw.occupancy[hashed_where[whos_hashed[i]]] = 0;
            }
            for (int i_mono = pivot_point+1; i_mono < n_mono; i_mono++){
                for (int j = 0; j < 3; j++){
                    trial_coord[i_mono][j] = coord[i_mono][j];
                }
            } 
            for (int i_mono = 0; i_mono < n_mono; i_mono++){
                add_to_hash_table(i_mono, trial_coord[i_mono]);
            }
        }

        

        

        //if(is_still_saw){
        for (int i_local = 0; i_local < n_mono/4; i_local++){
            /*if(uni_R(mt) < 0.5){
                single_bead_flip();
            }
            else{
                crankshaft_180();
            }*/
            int move = local_move_rand(mt);
            switch (move){
                case 0:
                    single_bead_flip();
                    break;
                case 1:
                    crankshaft_180();
                    break;
                case 2:
                    crankshaft_90_270(1);   // 90 degrees crankshaft
                    break;
                case 3:
                    crankshaft_90_270(-1);  // 270 degrees crankshaft
                    break;
                default:
                    std::cout << "INVALID LOCAL MOVE !!! \n\n\n";
            }
            
        }
        if (is_still_saw){n_pivots++;}
        compute_neighbour_list(trial_coord, trial_neighbours);
        trial_energy = compute_new_energy(trial_neighbours);
        //std::cout << trial_energy << "\n";
        delta_energy = trial_energy - energy;
        if (delta_energy <= 0){ acceptance = 1; }
        else {acceptance = exp(-delta_energy); }
        if (acceptance > uni_R(mt)){
            for (int i_mono = 0; i_mono < n_mono; i_mono++){
                for (int j = 0; j < 3; j++){
                    coord[i_mono][j] = trial_coord[i_mono][j];
                }
            }
            for (int i_mono = 0; i_mono < n_mono; i_mono ++){
                for (int j = 0; j < 7; j++){
                    neighbours[i_mono][j] = trial_neighbours[i_mono][j];
                }
            }
            energy = trial_energy;
            n_acc++;
        }
        //}
        // the next 3 lines clean up the hash table to make it ready again for use
        for (int i = 0; i < n_mono; i++){
            hash_saw.occupancy[hashed_where[i]] = 0;
        }
        
        spins_MC();  // here I run MC on the spin d.o.f's with the current polymer configuration
    }

    std::cout << "Number of successful pivot moves: "<< n_pivots << "\n";
    std::cout << "Number of actually accepted pivot moves: "<< n_acc << "\n";
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;
    myfile.open ("final_config_" + to_string(n_mono) + "_" + to_string(spin_coupling) + ".txt");
    myfile2.open ("e2e_dist_" + to_string(n_mono) + "_" + to_string(spin_coupling) + ".txt");
    myfile3.open ("energies_" + to_string(n_mono) + "_" + to_string(spin_coupling) + ".txt");
    myfile4.open ("final_spinconf_" + to_string(n_mono) + "_" + to_string(spin_coupling) +".txt");
    myfile5.open ("magnetization_" + to_string(n_mono) + "_" + to_string(spin_coupling) + ".txt");
    for (int i = 0; i < n_mono; i++){
        for(int j = 0; j < 3; j++){
            myfile << coord[i][j] << ' ';
        }
        myfile << "\n";
        myfile4 << spins[i] << "\n";
    }
    for (int i = 0; i < n_steps; i++){
        myfile2 << Ree2[i] << "\n";
        myfile3 << energies[i] << "\n";
        myfile5 << magnetization[i] << "\n";
    }
    myfile.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    myfile5.close();
}


/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/


int main(int argc, char* argv[]){
//int main(){
    stringstream string1(argv[1]);
    stringstream string2(argv[2]);
    stringstream string3(argv[3]);
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
    return 0;
}