# Likelihood-free inference for magnetic self-avoiding walks

## Current state

As of today (09/08/2021) the code written in "saw_with_ising.cpp" performs Monte Carlo simulations of a self avoiding walk, in which the monomers are also the spins of a Potts model with spatial nearest neighbours interactions.
The MC steps are perfomed as described in the following:
- Propose a global conformational move i.e. a pivot move (see the pivot algorithm from Madras & Sokal).
- If the proposed polymer is still self-avoiding, attemp an order N_monomers number of local moves (i.e. single bead flips and 2-beads crankshafts). Otherwise attemp the moves on the previous polymer configuration. 
- Accept or reject the resulting polymer configuration following the Metropolis-Hastings rule.
- Attempt an order N_monomers number of single spin flip MC moves.


## To Do list
