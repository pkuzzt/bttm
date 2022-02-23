#ifndef _MATRIX
#define _MATRIX
#include <stdio.h>

#define ERROR_INDEX_OVERFLOW 1
class coo_matrix {
    public:
        int *row_index;
        int *col_index;
        double *vals;
        int nnz;
        int nrow;
        int ncol;
        coo_matrix(int *row_index, int *col_index, double *vals, int nnz);
        ~coo_matrix();
};

class bcoo_matrix {
    public:
        int brow;
        int bcol;
        int nrow;
        int ncol;
        coo_matrix *subm;
        bcoo_matrix(int brow, int bcol, coo_matrix *A);
        ~bcoo_matrix();
};

void spmv(coo_matrix *A, double *x, double *b);
void spmv(bcoo_matrix *A, double *x, double *b);
#endif