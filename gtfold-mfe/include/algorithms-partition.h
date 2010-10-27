
#ifndef _ALGORITHMS_PARTITION_H
#define _ALGORITHMS_PARTITION_H

#ifdef __cplusplus
extern "C" {
#endif


void fill_partition_fn_arrays(int len, double** QB, double** Q, double** QM);
double **mallocTwoD(int r, int c);
void freeTwoD(double** arr, int r, int c);


#ifdef __cplusplus
}
#endif

#endif
