
#ifndef _ALGORITHMS_PARTITION_H
#define _ALGORITHMS_PARTITION_H

#ifdef __cplusplus
extern "C" {
#endif


void fill_partition_fn_arrays(int len, double** QB, double** Q, double** QM);
void printBasePairProbabilities(int n, int *structure, double **Q, double **QB);

double probabilityUnpaired(int index, int n, double **Q, double **QB);
double probabilityPaired(int i, int j, int n, double **Q, double **QB);

double **mallocTwoD(int r, int c);
void freeTwoD(double** arr, int r, int c);


#ifdef __cplusplus
}
#endif

#endif
