#include "matrix.hpp"
#include <algorithm>
#include <stdlib.h>
#include <omp.h>

using namespace std;
coo_matrix::coo_matrix(int *row_index, int *col_index, double *vals, int nnz) {
    this->row_index = (int*) malloc(nnz * sizeof(int));
    this->col_index = (int*) malloc(nnz * sizeof(int));
    this->vals = (double*) malloc(nnz * sizeof(double));
    this->nrow = 0;
    this->ncol = 0;
    for (int i = 0; i < nnz; i++) {
        this->row_index[i] = row_index[i];
        this->col_index[i] = col_index[i];
        this->vals[i] = vals[i];
        this->nnz = nnz;
        this->nrow = max(this->nrow, row_index[i]);
        this->ncol = max(this->ncol, col_index[i]);
    }
}

coo_matrix::~coo_matrix() {
    free(this->row_index);
    free(this->col_index);
    free(this->vals);
}

bcoo_matrix::bcoo_matrix(int brow, int bcol, coo_matrix *A) {
    this->brow = brow;
    this->bcol = bcol;    
    this->nrow = A->nrow;
    this->ncol = A->ncol;
    if (brow > A->nrow || bcol > A->ncol) {
        exit(ERROR_INDEX_OVERFLOW);
    }
    
    int nblock = brow * bcol;
    this->subm = (coo_matrix*) malloc(sizeof(coo_matrix) * nblock);
    for (int i = 0; i < nblock; i++) {
        this->subm[i].nnz = 0;
    }
    
    for (int i = 0; i < A->nnz; i++) {
        int irow = A->row_index[i] * brow / (A->nrow + 1);
        int icol = A->col_index[i] * bcol / (A->ncol + 1);
        int block_index = irow * bcol + icol;
        this->subm[block_index].nnz++;
    }

    int *index_list = (int*) malloc(sizeof(int) * nblock);
    for (int i = 0; i < nblock; i++) {
        index_list[i] = 0;
        int sub_nnz = this->subm[i].nnz;
        this->subm[i].row_index = (int*) malloc(sizeof(int) * sub_nnz);
        this->subm[i].col_index = (int*) malloc(sizeof(int) * sub_nnz);
        this->subm[i].vals = (double*) malloc(sizeof(double) * sub_nnz);
    }

    for (int i = 0; i < A->nnz; i++) {
        int irow = A->row_index[i] * brow / (A->nrow + 1);
        int icol = A->col_index[i] * bcol / (A->ncol + 1);
        int block_index = irow * bcol + icol;
        int local_index = index_list[block_index]++;
        coo_matrix* block_matrix = this->subm + block_index;
        block_matrix->row_index[local_index] = A->row_index[i];
        block_matrix->col_index[local_index] = A->col_index[i];
        block_matrix->vals[local_index] = A->vals[i];
    }
}

bcoo_matrix::~bcoo_matrix() {
    for (int i = 0; i < this->bcol * this->brow; i++) {
        free(this->subm[i].col_index);
        free(this->subm[i].row_index);
        free(this->subm[i].vals);
    }
}

void spmv(coo_matrix *A, double *x, double *b) {
    for (int i = 0; i < A->nrow; i++) {
        b[i] = 0;
    }
    for (int i = 0; i < A->nnz; i++) {
        int i1 = A->row_index[i];
        int i2 = A->col_index[i];
        double val = A->vals[i];
        b[i1] += val * x[i2];
    }
}

void spmv(bcoo_matrix *A, double *x, double *b) {
    for (int i = 0; i < A->nrow; i++) {
        b[i] = 0;
    }
    for (int block = 0; block < A->bcol * A->brow; block++) {
        coo_matrix *M = A->subm + block;
        for (int i = 0; i < M->nnz; i++) {
            int i1 = M->row_index[i];
            int i2 = M->col_index[i];
            double val = M->vals[i];
            b[i1] += val * x[i2];
        }
    }
}