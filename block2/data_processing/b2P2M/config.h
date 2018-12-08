


#ifndef F4_config_h
#define F4_config_h






// ********************* inlcudes ************************

// standard lib
#include <time.h>
#include <limits.h>
#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
//#include <sys/malloc.h>

// pthread
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



// m4ri
#include <m4ri/config.h>
#include <m4ri/m4ri.h>

// cudd
#include <util.h>
#include <cuddInt.h>
#include <cudd.h>







//! maximal degree
#define DEGREE_LIMIT 8 //! = 0 mod 8




// ! type redefinition
typedef unsigned long gb_word_t;
typedef unsigned short var_subscript_t; //! subscripts of xi. 0 <= i < 2^16
typedef unsigned char data_t; //! data in sparse matrices.
typedef char deg_t; //! degree
typedef int idrc_t; //! size of rows and cols in a matrix, or size of monomials.
typedef short mindex_t; //! index of matrix.
typedef int blockindex_t; //! index of block
typedef short blockpos_t; //! position in a block.
typedef char flag_t; //! flags or status
typedef long count_t;
//typedef unsigned long hashvalue_t;



//! the maximal addmissible number of variables is MON_LEN * WORD_SIZE
const unsigned char MON_LEN = 2;
const var_subscript_t WORD_SIZE = sizeof(gb_word_t) * 8;
const idrc_t WORD_MASK = WORD_SIZE - 1;
const unsigned char WORD_OFFSET  = 6;
const var_subscript_t VAR_MASK = WORD_SIZE - 1;
const unsigned char VAR_OFFSET  = 6;


const var_subscript_t VAR_LIMIT = MON_LEN * WORD_SIZE;
typedef gb_word_t fullmon_t[MON_LEN];
typedef var_subscript_t mon_t[DEGREE_LIMIT];

const mindex_t MATRIX_LIMIT  = 32;
const idrc_t DYNPOSTOID_LIMIT = MATRIX_LIMIT * VAR_LIMIT;






//! useful functions
#define GB_MAX(a, b) (a > b ? a : b)
#define GB_MIN(a, b) (a < b ? a : b)


#define INFO_PR(s) {printf("\n************ %s *************\n\n", s);}


#define ERROR(s) {printf("\n ERROR !!!!!!------------- %s ----------------\n\n", s); exit(1);}

#define ONEAPPEAR {printf("\n+++++++ 1 appears +++++++\n\n"); return;}

#define PR_ENTER {printf("\n");}


#define MAT_PR(pivot, nonpivot, spoly) {printf("~~~~~~~~~~ matrix (%d x %d) information [A: %d x %d] [B: %d x %d] [C: %d x %d] [ D: %d x %d]  ~~~~~~~~~~\n", (pivot) + (spoly), (pivot) + (nonpivot), (pivot), (pivot), (pivot), (nonpivot), (spoly), (pivot), (spoly), (nonpivot));}

#define LOOPINFO_PR(s, d) {char buf[80]; sprintf(buf, s, d); printf("\n***************** %s ***************** \n\n", buf);}


#if DEBUG_MODE

#define LOOP_PR(s) {printf("~~~~~~~~~~ %s ~~~~~~~~~~\n", s);}
#define LOOPIN_PR(s, d) {char buf[80]; sprintf(buf, s, d); printf("~~~~~~~~~~ %s ~~~~~~~~~~\n", buf);}
#else
#define LOOP_PR(s) {}
#define LOOPIN_PR(s, d) {}

#endif






















#endif
