
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "algorithms.h"
#include "algorithms-partition.h"
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

    // Boltzmann constant (R) * Standard 37C temperature (T in Kelvin)
    double RT = 0.00198721 * 310.15;

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

/**
 * @param n Length of the RNA strand
 * @param structure Array of what pairs with what.  structure[i] = j means
 *                  that the nucleotide at index i is paired with that at
 *                  index j. structure[i] = 0 means the nucleotide is
 *                  unpaired.
 */
void fillBasePairProbabilities(int length, int *structure, double **Q, double **QB, double **QM, double**P) {
	int d, l, h, i, j;
	double tempBuffer;
    	// multiConst[3] is a global variable with 3 values: a, b, c for the
	// experimental constants
 	int a = multConst[0]; // a is an offset penalty for multiloops
	int b = multConst[2]; // b is penalty for multiloop branches, one per branchesch
	int c = multConst[1]; // Penalty for single stranded nucleotides in thee multiloops
	for(d = length; d>3; d--){
		for(h = 1; h + d <= length; h++){
			for(l = h + d; l<=length; l++){
				P[h][l] = (QB[h][l] / Q[1][length]);  //Prob. h,l is an exterior base pair
				if(h > 1)
					P[h][l] *= Q[1][h-1];
				if(l < length)
					P[h][l] *= Q[l+1][length];
				for(i = 1; i < h; i++)
					for(j = l+1; l <= length; l++){ //Now For internal loops
						tempBuffer = P[i][j]*QB[h][l]/QB[i][j];
						if(i == h-1 && j == l+1) //of which stacked pairs are a special case
							tempBuffer *= exp(-eS(i,j,h,l)/RT);
						else
							tempBuffer *= exp(-eH(i,j,h,l)/RT);

						P[i][j] += tempBuffer;

						tempBuffer = 0; // Start over for multiloops
						if(j - l > 3)
							tempBuffer += exp(-((h-i-1)*c/RT * QM[l+1][j-1]));
						if(h - i > 3)
							tempBuffer += exp(-((j-l-1)*c/RT) * QM[i+1][h-1]);
						if(j - l > 3 && h -i > 3)
							tempBuffer += QM[i+1][h-1] * QM[l+1][j-1];

						tempBuffer *= P[i][j] * QB[h][l] / QB[i][j] * exp(-(a+b)/RT);

						P[i][j] += tempBuffer;
					}
			}
		}
	}


}

/**
 * Implement this expression:
 *
 *      n
 * 1 - sum (Q_{1,i-1} * QB_{i,j} * Q_{j+1,n})
 *      1
 * ------------------------------------------
 *             Q_{1,n}
 *
 * @param index The nucleotide to consider
 * @param n The length of the RNA strand
 * @param Q The Q matrix
 * @param QB The QB matrix
 */
double probabilityUnpaired(int index, int n, double **Q, double **QB) {

    double p = 0;

    int i;
    for(i=1; i<=n; ++i)
        if(i != index)
            p += probabilityPaired(index, i, n, Q, QB);

    return 1 - p;
}

/**
 *
 * Q_{1,i-1} * QB_{i,j} * Q_{j+1,n}
 * --------------------------------
 *            Q_{1,n}
 *
 * @param i The first nucleotide to consider
 * @param j The second nucleotide to consider
 * @param n The length of the RNA strand
 * @param Q The Q matrix
 * @param QB The QB matrix
 */
double probabilityPaired(int i, int j, int n, double **Q, double **QB) {

    // ensure that the invariant i<=j holds
    int newi = MIN(i,j);
    int newj = MAX(i,j);
    i = newi;
    j = newj;

    // the first and the last nucleotides pair
    if(i == 1 && j == n)
        return QB[i][j] / Q[1][n];

    // j pairs with the first nucleotide
    else if(i == 1)
        return QB[i][j] * Q[j+1][n] / Q[1][n];

    // i pairs with last nucleotide
    else if(j == n)
        return Q[1][i-1] * QB[i][j] / Q[1][n];

    // neither pairing isn't with either first or last nucleotide
    else
        return Q[1][i-1] * QB[i][j] * Q[j+1][n] / Q[1][n];
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

