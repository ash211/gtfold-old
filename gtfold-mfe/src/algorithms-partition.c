
#include <math.h>

#include "algorithms.h"
#include "data.h"


// double[][] QB;
// double[][] Q;
// double[][] QM;


// Based on pseudocode in figure 6 in 
//
// A Partition Function Algorithm for Nucleic Acid Secondary Structure
// Including Pseudoknots
// Dirks&Pierce 2003
void fillArrays(int len, double** QB, double** Q, double** QM) {

    // Boltzmann constant (R) * Standard 37C temperature (T)
    double RT = 8.31447215 * 310.15;

    // multiConst[3] is a global variable with 3 values: a, b, c for the
    // experimental constants
    int a = multConst[0]; // a is an offset penalty for multiloops
    int b = multConst[2]; // b is penalty for multiloop branches, one per branchesch
    int c = multConst[1]; // Penalty for single stranded nucleotides in thee multiloops


    // initialize arrays
    for(int i=0; i<len; ++i)
        for(int j=0; j<len; ++j)
            QB[i][j] = Q[i][j] = QM[i][j] = 0;

    for(int i=1; i<len; ++i)
        Q[i][i-1] = 1;


    // fill in values in the array
    for(int l=1; i<=N; ++i) {
        for(int i=1; i<= N-l+1; ++i) {
            int j = i+l-1;

            // QB recursion
            // NOTE: eH returns an integer encoded as fixed point.  So a
            // return value of 115 represents raw value 115/100 = 1.15
            QB[i][j] = exp(-eH(i,j)/100.0/RT);

            for(int d=i+1; d<=j-4; ++d) {
                for(int e=d+4; e<=j-1; ++e) {
                    QB[i][j] += exp(-eL(i,j,d,e)/100.0/RT)*QB[d][e];

                    QB[i][j] += QM[i+1][d-1]*QB[d][e] *
                        exp(-(a + 2b + c*(j-e-1))/RT);
                }
            }

            // Q, QM recursions
            Q[i][j] = 1;
            for(int d=i; d<=j-4; ++d) {
                for(int e=d+4; e<=j; ++e) {
                    Q[i][j] += Q[i][d-1]*QB[d][e];
                    QM[i][j] += exp(-(b+c*(d-i)+c*(j-e))) * QB[d][e];
                    QM[i][j] += QM[i][d-1] * QB[d][e] * exp(-(b+c*(j-e))/RT);
                }
            }
        }
    }
}
