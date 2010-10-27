
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
/**
 *
 * @param N The length of the RNA strand
 * @param QB Matrix
 * @param Q Matrix
 * @param QM Matrix
 *
 */
void fill_partition_fn_arrays(int len, double** QB, double** Q, double** QM) {

    // Boltzmann constant (R) * Standard 37C temperature (T)
    double RT = 8.31447215 * 310.15;

    // multiConst[3] is a global variable with 3 values: a, b, c for the
    // experimental constants
    int a = multConst[0]; // a is an offset penalty for multiloops
    int b = multConst[2]; // b is penalty for multiloop branches, one per branchesch
    int c = multConst[1]; // Penalty for single stranded nucleotides in thee multiloops

    // loop iterators
    int i,j,d,e,l;

    fprintf(stdout, "Initializing PF arrays. . .");
    // initialize arrays
    for(i=0; i<len; ++i) {
        for(j=0; j<len; ++j) {
            QB[i][j] = 0;
            Q[i][j] = 0;
            QM[i][j] = 0;
        }
    }


    for(i=1; i<len; ++i)
        Q[i][i-1] = 1;

    fprintf(stdout, " done.\n");
    fprintf(stdout, "Filling PF arrays. . .");

    // fill in values in the array
    for(l=1; l<=len; ++l) {
        fprintf(stdout, "Running with l=%d\n", l);
        for(i=1; i<= len-l+1; ++i) {

            int j = i+l-1;

            // QB recursion
            // NOTE: eH returns an integer encoded as fixed point.  So a
            // return value of 115 represents raw value 115/100 = 1.15
            QB[i][j] = exp(-eH(i,j)/100.0/RT);

            for(d=i+1; d<=j-4; ++d) {
                for(e=d+4; e<=j-1; ++e) {
                    QB[i][j] += exp(-eL(i,j,d,e)/100.0/RT)*QB[d][e];

                    QB[i][j] += QM[i+1][d-1]*QB[d][e] *
                        exp(-(a + b + c*(j-e-1))/RT);
                }
            }

            // Q, QM recursions
            Q[i][j] = 1;
            for(d=i; d<=j-4; ++d) {
                for(e=d+4; e<=j; ++e) {
                    Q[i][j] += Q[i][d-1]*QB[d][e];
                    QM[i][j] += exp(-(b+c*(d-i)+c*(j-e))) * QB[d][e];
                    QM[i][j] += QM[i][d-1] * QB[d][e] * exp(-(b+c*(j-e))/RT);
                }
            }
        }
    }
}


double **mallocTwoD(int r, int c) {
    double** arr = (double **)malloc(r*sizeof(double));
    int i;
    for(i=0; i<r; i++) {
        arr[i] = (double *)malloc(c*sizeof(double));

        // failed allocating a row, so free all previous rows, free the main
        // array, and return NULL
        if(arr[i] == NULL) {
            int j;
            for(j=0; j<i; j++)
                free(arr[j]);
            free(arr);
            return NULL;
        }
    }

    return arr;
}

void freeTwoD(double** arr, int r, int c) {
    int i;
    for(i=0; i<r; i++)
        free(arr[i]);

    free(arr);
}

