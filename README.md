# Likelihood-free inference for magnetic self-avoiding walks

## Current state

As of today (09/08/2021) the code written in "saw_with_ising.cpp" performs Monte Carlo simulations of a self avoiding walk, in which the monomers are also the spins of a Potts model with spatial nearest neighbours interactions.
The MC steps are perfomed as described in the following:
- Propose a global conformational move i.e. a pivot move (see the pivot algorithm from Madras & Sokal).
- If the proposed polymer is still self-avoiding, attemp an order N_monomers number of local moves (i.e. single bead flips and 2-beads crankshafts). Otherwise attemp the moves on the previous polymer configuration. 
- Accept or reject the resulting polymer configuration following the Metropolis-Hastings rule.
- Attempt an order N_monomers number of single spin flip MC moves.

06/09/2021
- Now the code performs properly the MMC algorithm and a check has been made by comparing plots obtained in the zero field scenario with a paper by Garel et al (1999) about magnetic polymers
- Currently working on a very rudimental implementation of the Wood SL method on synthetic data where the fields are just an unknown weight that multiplies a single sequence feature, i.e. the CpG ratio of the fragment that is represented by a single monomer


## To Do list
- ✔ Simulations with 0/1 spins instead of a classical Ising model (ChIP-seq should give positive numbers)
- Download some ChIP-seq data and try to understand how the files are structured. Look up the peak-calling thing that Guido mentioned.
- Think about how to summarize the data.
- Write functions to compute summary statics.
- Write the code to perform the Monte Carlo in the space of the parameters, with a method such as synthetic likelihood (Wood, 2010).
