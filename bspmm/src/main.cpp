#include "matrix.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main(int argc, char* argv[]) {
    int nrow = 2000000;
    int ncol = 2000000;
    int nnz = 100000000;
    int* row_index = (int*) malloc(sizeof(int) * nnz);
    int* col_index = (int*) malloc(sizeof(int) * nnz);
    double *x = (double*) malloc(sizeof(double) * ncol);
    double *b = (double*) malloc(sizeof(double) * nrow);
    double *vals = (double*) malloc(sizeof(double) * nnz);
    for (int i = 0; i < ncol; i++) {
        x[i] = drand48();
    }
    for (int i = 0; i < nnz; i++) {
        row_index[i] = rand() % nrow;
        col_index[i] = rand() % ncol;
        vals[i] = drand48();
    }
    double t1, t2;
    coo_matrix A(row_index, col_index, vals, nnz);
    bcoo_matrix B(100, 100, &A);
    int times = 10;
    t1 = omp_get_wtime();
    for (int t = 0; t < times; t++) {
        spmv(&A, x, b);
    }
    t2 = omp_get_wtime();
    fprintf(stdout, "coo format: %f\n", t2 - t1);

    t1 = omp_get_wtime();
    for (int t = 0; t < times; t++) {
        spmv(&B, x, b);
    }
    t2 = omp_get_wtime();
    fprintf(stdout, "blocked coo format: %f\n", t2 - t1);
    return 0;
}
