# Likelihood-free inference for magnetic self-avoiding walks

## Current state

As of today (09/08/2021) the code written in "saw_with_ising.cpp" performs Monte Carlo simulations of a self avoiding walk, in which the monomers are also the spins of a Potts model with spatial nearest neighbours interactions.
The MC steps are perfomed as described in the following:
- Propose a global conformational move i.e. a pivot move (see the pivot algorithm from Madras & Sokal).
- If the proposed polymer is still self-avoiding, attemp an order N_monomers number of local moves (i.e. single bead flips and 2-beads crankshafts). Otherwise attemp the moves on the previous polymer configuration. 
- Accept or reject the resulting polymer configuration following the Metropolis-Hastings rule.
- Attempt an order N_monomers number of single spin flip MC moves.


## To Do list
### Part 1: simulator
- ✔ Put everything you have written until now in a library (magpol.a) and a main program that calls this library (otherwise code becomes too messy and impossible to use going forward) 
- Instead of taking simulation parameters from the command line when the code is launched, write the code in such a way that all the required info comes from an input text file.
- Separate the part that writes the outputs on files from the "saw_MC::run()" method.
- Instead of initializing the polymer as a straight rod and the spin configurations randomly, you should also allow the possibility to read a configuration from a file. This is required first of all because when you do the likelihood-free inference part, you don't want to wait for a long equilibration for each set of model parameters, but instead you would want to start from a configuration that is already close to equilibration. The second reason is that as it is now the MC "dynamics" is metastable. Especially in the case of zero external fields and strong coupling between the spins almost all global polymer moves are rejected. The local moves are accepted more easily but they do not change the end-to-end distance.
- Implement the parallel tempering/multiple markov chain method to try to overcome this metastable behaviour. (many questions: how do you choose the temperature range? How many temperatures? How often to attempt a swap? ecc...)
### Part 2: inference and data
- Download some ChIP-seq data and try to understand how the files are structured. Look up the peak-calling thing that Guido mentioned.
- Think about how to summarize the data.
- Write the code to perform the Monte Carlo in the space of the parameters, with a method such as synthetic likelihood (Wood, 2010).
