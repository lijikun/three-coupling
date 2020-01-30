Numerically calculates the transfer Hamiltonian matrix T of a particle with spin Se hopping between two identical sites with angular momentum Ss. Can be modified to calculate more general spin Hamiltonian matrices.

Depends on GNU Scientific Library (GSL) and optionally OpenMP.

To compile:

    g++ --std=c++11 -O2 -fopenmp -lgsl -lgslcblas -lm couple-three.cc -o couple-three

Usage: 

    ./couple-three 2Ss 2Se 2S1 [2S2 2S3 ...]

Where 2Ss is twice the angular momentum of a site, and 2Se is twice the spin of the hopping particle, while 2S1, 2S2, 2S3, ... are twice total spins for the subspaces to be calculated -- each S subspace has 2S+1 basis wavefunctions, i.e. the states |S, M> where M = -S, -S+1, ..., +S-1, +S.

Because all angular momentum quantum numbers are integers or half integers, the program takes twice their values (2S) as command line parameters, which are expected to be integers. For such a system the possible value of total 2S ranges from 0 to 2(2Ss) + 2Se and must have the same parity (i.e. being both odd or even) as 2Se. Invalid parameter value(s) will cause an error message, while duplicate input(s) will be silently ignored.

One can use [GNU seq](https://www.gnu.org/software/coreutils/manual/html_node/seq-invocation.html) to generate a range of possible values, e.g.: `$(seq 1 2 31)` will be interpreted as sequence of `1 3 5 ... 31` in the command line.
