// Calculates the hopping Hamiltonian matrix of a system where a particle "e" with angular momentum Se 
// hops between two identical sites "a" and "b", each site has angular momentum Ss.
// S is the total angular momentum with many possible values, but only calculate the ones provided by the user.
// All angular momentum are integers or half integers, therefore twice values (2Ss, 2Se, 2S, etc.) are stored as int.
//
// To compile: g++ --std=c++11 couple-three.cc -lgsl -lgslcblas -lm -fopenmp -o couple-three

#include <cstdio>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

double cg_coef(int, int, int, int, int, int);

int main(int argc, char* argv[]) {
    std::vector<int> twoS_values;
    std::vector<int> mat_index_bases;
    std::map<int,int> twoS_index;
    int twoSs, twoSe, temp_2S;
    // Handles input. Need at least 3 integers.
    // Format: ./couple-three 2Sa 2Se 2S1 2S2 2S3 ...
    // Todo: Allow input without S1, S2, etc. (calculate for all possible S by default)
    // Todo: Allow input format like Slow-Shigh, with intervals of 2. Need to check parity.
    try {
        if (argc < 4) throw 1;
        twoSs = std::stoi(argv[1]);
        twoSe = std::stoi(argv[2]);
        if (twoSs < 1 || twoSe < 1) throw 2;
        mat_index_bases.push_back(0);
        for (int i = 3; i < argc; ++i) {
            temp_2S = std::stoi(argv[i]);
            // 0 <= S <= 2Ss+Se and 2S should have the same parity with 2Se.
            if ( temp_2S > twoSs * 2 + twoSe || temp_2S < 0 || (temp_2S & 1) != (twoSe & 1) ) throw 3;
            // Uses a map to provide reverse lookup from an 2S value to its index in twoS_values.
            // If an S is listed more than once, inserts after the 1st time will be skipped silently.
            std::pair<std::map<int,int>::iterator, bool> try_insert;
            try_insert = twoS_index.insert({temp_2S, i-3});
            if (try_insert.second) {
                twoS_values.push_back(temp_2S);
                // Each S adds 2S+1 states.
                mat_index_bases.push_back(mat_index_bases.back() + temp_2S + 1);
            }
        }
    }
    catch (...) {
        printf("Erroneous parameters found. Correct usage: ./couple-three 2Ss 2Se 2S1 [2S2 2S3 ...]\n\n");
        printf("At least 3 non-negative integer parameters needed satisfying the rules of angular momentum coupling.\n");
        printf("Ss is the angular momentum for the two identical sites and Se is that of the hopping particle.\n");
        printf("S1, S2, etc specify total angular momentum subspaces to be calculated.\n");
        printf("To specify a wide range, use $(seq start step end), e.g.: ./couple-three 15 1 $seq(1 2 31)\n");
        return 1;
    }
    // The last element of mat_index_bases is the sum of all 2S+1, hence the dimensionality of the matrix.
    std::vector<std::vector<double>> t_matrix(mat_index_bases.back(), std::vector<double>(mat_index_bases.back(), 0.0));
    
    // Calculates each matrix element.
    //printf("2Ma, 2Mb, 2Me, 2Sinta, 2Sintb, 2Sa, 2Sb, i, j, mat_elem\n");
    int twoMa, twoMb, twoMe;
    #pragma omp parallel for collapse(3) private(twoMa, twoMb, twoMe) shared(t_matrix)
    for (twoMa = -twoSs; twoMa <= twoSs; twoMa +=2) {
        for (twoMb = -twoSs; twoMb <= twoSs; twoMb +=2) {
            for (twoMe = -twoSe; twoMe <= twoSe; twoMe +=2) {
                //#ifdef _OPENMP
                //if (twoMa == twoSs && twoMb == twoSs & twoMe == twoSe)
                //    printf("Ran with %d threads.\n", omp_get_num_threads());
                //#endif
                // Total M.
                int twoM = twoMa + twoMb + twoMe;
                // Sa & Se couple into |Sinta, Ma+Me> state
                // Sb & Se couple into |Sintb, Mb+Me> state
                for (int twoSinta = std::max(abs(twoMa+twoMe), abs(twoSs-twoSe)); twoSinta <= twoSs+twoSe; twoSinta +=2) {
                    for (int twoSintb = std::max(abs(twoMb+twoMe), abs(twoSs-twoSe)); twoSintb <= twoSs+twoSe; twoSintb +=2) {
                        // Iterates over possible S values to see if the two total states, |Sa, M> and |Sb, M>,
                        // are possible. If yes, calculates the contribution to matrix element <Sa, M| T |Sb, M>
                        for (auto&& twoSa : twoS_values ) {
                            for (auto&& twoSb : twoS_values) {
                                if (twoSa >= abs(twoM) && twoSa >= abs(twoSinta-twoSs) && twoSa <= twoSinta+twoSs
                                    && twoSb >= abs(twoM) && twoSb >= abs(twoSintb-twoSs) && twoSb <= twoSintb+twoSs) {
                                    int i = mat_index_bases[twoS_index[twoSa]] + (twoM+twoSa)/2;
                                    int j = mat_index_bases[twoS_index[twoSb]] + (twoM+twoSb)/2;
                                    double mat_elem_temp = cg_coef(twoSinta, twoMa+twoMe, twoSs, twoMb, twoSa, twoM) 
                                        * cg_coef(twoSs, twoMa, twoSe, twoMe, twoSinta, twoMa+twoMe)
                                        * cg_coef(twoSintb, twoMb+twoMe, twoSs, twoMa, twoSb, twoM)
                                        * cg_coef(twoSs, twoMb, twoSe, twoMe, twoSintb, twoMb+twoMe);
                                    #pragma omp atomic
                                    t_matrix[i][j] += mat_elem_temp;
                                    //printf("%d %d %d %d %d %d %d %d %d %f\n", twoMa, twoMb, twoMe, twoSinta, twoSintb, twoSa, twoSb, i, j, mat_elem_temp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Output
    for (int i = 0; i < mat_index_bases.back(); ++i) {
        for (int j = 0; j < mat_index_bases.back(); ++j) 
            printf("%9f ", t_matrix[i][j]);
        printf("\n");
    }

    return 0;
}

// Clebsch-Gordan coefficients, wrapper for the GSL 3j-symbol function.
double cg_coef(int two_j1, int two_m1, int two_j2, int two_m2, int two_j3, int two_m3) {
    return ( (1&(two_j2-two_j1-two_m3) / 2) ? -1.0 : 1.0 ) * sqrt(double(two_j3) + 1.0) 
        * gsl_sf_coupling_3j( two_j1, two_j2, two_j3, two_m1, two_m2, -two_m3 );
}
