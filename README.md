Numerically calculates the transfer Hamiltonian matrix T of a particle with spin Se hopping between two identical sites with angular momentum Ss. Can be modified to calculate more general spin Hamiltonian matrices.

Depends on GNU Scientific Library (GSL) and optionally OpenMP.

To compile:

    g++ --std=c++11 -O2 -fopenmp -lgsl -lgslcblas -lm couple-three.cc -o couple-three

Usage: 

    ./couple-three 2Ss 2Se 2S1 [2S2 2S3 ...]

Because all angular momentum quantum numbers (S) are integers or half integers, the program takes twice their values (2S) as command line parameters. For such a system possible values total S ranges from 0 to 2Ss + Se and has the same parity as Se. Invalid parameter values will cause an error message, while duplicate input will be silently ignored.

One can use [GNU seq](https://www.gnu.org/software/coreutils/manual/html_node/seq-invocation.html) to generate a range of possible values, e.g.: `$(seq 1 2 31)` will be interpreted as sequence of `1 3 5 ... 31` in the command line.
