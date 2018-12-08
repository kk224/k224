
#include "evaluation.h"

void evaluation_init(evaluation_t *eval,zddpoly_t *zddp){
	    
    memset(eval->base_mat,0,sizeof(eval->base_mat));
    memset(eval->constant_array,0,sizeof(eval->base_mat));
    memset(eval->var_array,0,sizeof(eval->var_array));
    memset(eval->quadratic_mat,0,sizeof(eval->quadratic_mat));
    memset(eval->nonlinear_mat,0,sizeof(eval->nonlinear_mat));
    for (idrc_t i = 0; i < NONLINEAR_NUM; i++) {
        zdd_to_matrix_total(eval, i, zddp->zdd, zddp->polys[i], zddp->old_to_new);

    }
    
    for (idrc_t i = NONLINEAR_NUM; i < zddp->size_of_polys; i++) {

        zdd_to_column(eval, i-NONLINEAR_NUM, zddp->zdd, zddp->polys[i], zddp->old_to_new);
    }
    
//    base_mat_print(eval);
}


void zdd_to_matrix_total(evaluation_t *eval,int mid,DdManager *zdd, DdNode *p, mzp_t *old_to_new){
    
    DdNode *t, *e;
    DdNode *one = DD_ONE(zdd);
    DdNode *zero = DD_ZERO(zdd);
    int var_all = PARAM_NUM+VAR_NUM;
    int row = 0;
    for (int i = 0; i< PARAM_NUM + VAR_NUM + 1; i++) {
        
        if (p == zero){
            continue;
        }else if (p == one){
            
            eval->nonlinear_mat[mid][var_all][var_all] = 1;
        }else{
            row = old_to_new->values[p->index];
            
            t = cuddT(p);
            if(t == one){
               
                //todo var_array
                eval->nonlinear_mat[mid][row][var_all] = 1;
                
                
            }else if (t == zero){
                printf("\nZDD is WRONG!!!\n\n\n");
            }else{
                eval->nonlinear_mat[mid][row][old_to_new->values[t->index]] = 1;
                
                e = cuddE(t);
                
                while (e != zero && e != one) {
                    
                    eval->nonlinear_mat[mid][row][old_to_new->values[e->index]] =1;
                    
                    e = cuddE(e);
                }
                
                if (e == one) {
                    eval->nonlinear_mat[mid][row][var_all] = 1;
                }
            }
            p = cuddE(p);
        }
        
    }
    

}

void zdd_to_column(evaluation_t *eval,int col, DdManager *zdd, DdNode *p, mzp_t *old_to_new){
    
    DdNode *t, *e;
    DdNode *one = DD_ONE(zdd);
    DdNode *zero = DD_ZERO(zdd);
    
    
    unsigned long base = 0;
    unsigned long constant = 0;
    int row = 0;
    for (int i = 0; i< PARAM_NUM + VAR_NUM + 1; i++) {
        base = 0;
        if (p == zero){
            continue;
        }else if (p == one){
            
            eval->quadratic_mat[col][PARAM_NUM][PARAM_NUM] = 1;
        }else{
            
            if (p->index< 100) {  //Ci*(Ci's+Xi's)
                
                row = old_to_new->values[p->index];
                
                t = cuddT(p);
                if(t!= zero && t!= one ){
                    if (t->index > 100) {   //Xi's
                        
                        base ^= ((unsigned long)1 << (PARAM_NUM + VAR_NUM - 1 - old_to_new->values[t->index]));
                        
                    }else{  //Ci's
                        
                        eval->quadratic_mat[col][row][old_to_new->values[t->index]] =1;
                    }
                    
                    e = cuddE(t);
                    while (e != zero && e != one) {
                                                      
                        if(e->index > 100){
                            
                            base ^= ((unsigned long)1<< (PARAM_NUM + VAR_NUM - 1 -old_to_new->values[e->index]));
                            
                            
                            
                        }else{
                                                      
                            eval->quadratic_mat[col][row][old_to_new->values[e->index]] =1;
                                                      
                        }
                        e = cuddE(e);
                    }
                    if (e == one) {
                        
                        //todo var_array
                        eval->quadratic_mat[col][row][PARAM_NUM] = 1;
                        
                        
                    }
                    
                }else if(t == one){
                    
                    //todo var_array
                    eval->quadratic_mat[col][row][PARAM_NUM] = 1;
                    
                    
                }else if (t == zero){
                    printf("\nZDD is WRONG!!!\n\n\n");
                }
                
                
                eval->base_mat[row][col] = base;
                
                
                
                
            }else{  //Xi*(Xi's)
                
                if(cuddT(p)!= one){
                    printf("not a linear polys after evaluation!!!!\n\n");
                    
                    break;
                }else{
                    
                    constant ^= ((unsigned long)1<< (PARAM_NUM + VAR_NUM - 1 -old_to_new->values[p->index]));
                }
                
                
                
                
                
            }
            p = cuddE(p);
            
        }
        eval->constant_array[col] = constant;
        
    }
    
}

void base_mat_print(evaluation_t *eval){
    
    printf("\n------------ base  matrix ---------------\n");
    
    for (int i = 0; i< PARAM_NUM; i++) {
        for (int j = 0; j < LINEAR_NUM; j++) {
            
            printf("%lu ",eval->base_mat[i][j]);
        }
        printf("\n");
    }
    
    printf("\n------------ constant array ---------------\n");
    
    for (int i = 0; i < LINEAR_NUM; i++) {
        printf("%lu\n",eval->constant_array[i] );
    }
    
    printf("\n------------ quadratic matrices ---------------\n");
    for (int i = 0; i < LINEAR_NUM; i++) {
        printf("%d polys\n",i);
        for (int j = 0; j < PARAM_NUM+1; j++) {
            for (int k = 0; k < PARAM_NUM+1; k++) {
                printf("%d ",eval->quadratic_mat[i][j][k]);
            }
            printf("\n");
        }
    }
    
    printf("\n------------ nonlinear matrices ---------------\n");
    for (int i = 0; i < NONLINEAR_NUM; i++) {
        printf("%d nonlinear polys\n",i);
        for (int j = 0; j < PARAM_NUM+ VAR_NUM+1; j++) {
            for (int k = 0; k < PARAM_NUM +VAR_NUM +1; k++) {
                printf("%d ",eval->nonlinear_mat[i][j][k]);
            }
            printf("\n");
        }
    }

    
}


