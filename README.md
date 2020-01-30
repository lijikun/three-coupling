Numerically calculates the transfer Hamiltonian matrix T. Can be modified to calculate more general spin Hamiltonian matrices.

Depends on `gsl` and optionally `openmp`.

To compile:

    g++ --std=c++11 -O2 -fopenmp -lgsl -lgslcblas -lm couple-three.cc -o couple-three

Usage: 

    ./couple-three 2Ss 2Se 2S1 [2S2 2S3 ...]

Because all angular momentum quantum numbers (S) are integers or half integers, the program takes twice values (2S) as integer inputs. For such a system possible values total S ranges from 0 to 2S_s + S_e and has the same parity as S_e. One can use [GNU seq](https://www.gnu.org/software/coreutils/manual/html_node/seq-invocation.html) to generate a range of possible values, i.e. `$(seq 1 2 31)` will generate a sequence of `1 3 5 ... 31` and append it to the command line.
