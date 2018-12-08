
#ifndef F4_evaluate_h
#define F4_evaluate_h

#include <stdio.h>
#include "zdd.h"

#define PARAM_NUM 37
#define VAR_NUM  187
#define LINEAR_NUM 190
#define NONLINEAR_NUM 0

typedef long value_t; // to save values of variables.

//typedef unsigned long constpart_t; // the part with no parameters.
typedef unsigned long linearpart_t; // to save 32 unknowns and 1 contants.
typedef unsigned long squarepart_t;

typedef unsigned long oripoly_t;

struct evaluation_t{
    unsigned long base_mat[PARAM_NUM][LINEAR_NUM][3];
    unsigned long constant_array[LINEAR_NUM][3];
    bool quadratic_mat[LINEAR_NUM][PARAM_NUM+1][PARAM_NUM+1];
    bool nonlinear_mat[NONLINEAR_NUM][PARAM_NUM+VAR_NUM+1][PARAM_NUM+VAR_NUM+1];
    
    
};


void evaluation_init(evaluation_t *eval,zddpoly_t *zddp);
void zdd_to_column(evaluation_t *eval,int col, DdManager *zdd, DdNode *p, mzp_t *old_to_new);
void base_mat_print(evaluation_t *eval);
void zdd_to_matrix_total(evaluation_t *eval,int mid,DdManager *zdd, DdNode *p, mzp_t *old_to_new);







/// functions

static inline void binary_print(value_t val, int len) {
    
    for (int i = 0; i < len; i++) {
        if (val & ((value_t)1 << i)) {
            printf("1");
        } else {
            printf("0");
        }
        
        if ((i + 1) % 5 == 0) {
            printf(" ");
        }
    }
    
}


static inline int largestpos(value_t val, int len) {
    
    CHECK_GT(val, 0)
    
    for (int i = len - 1; i >= 0; i--) {
        if (val & ((value_t)1 << i)) {
            return i;
        }
    }
    
    return -1;
}



static inline void print_status(linearpart_t *status, int poly_num, int word_len) {
    
    for (int i = 0; i < poly_num; i++) {
        binary_print(status[i], word_len);
        
        if (i < poly_num - 1) {
            printf("\n");
        } else {
            printf("\n");
        }
    }
    
}





#endif
