
#include "evaluate.h"



void evaluation_init(evaluation_t *eval,zddpoly_t *zddp){

    memset(eval->base_mat,0,sizeof(eval->base_mat));
    memset(eval->constant_array,0,sizeof(eval->base_mat));
    memset(eval->quadratic_mat,0,sizeof(eval->quadratic_mat));
    memset(eval->nonlinear_mat,0,sizeof(eval->nonlinear_mat));
    for (idrc_t i = 0; i < NONLINEAR_NUM; i++) {
        zdd_to_matrix_total(eval, i, zddp->zdd, zddp->polys[i], zddp->old_to_new);

    }
    
    for (int i = 0; i < LINEAR_NUM; i++) {
        for (int j = 0; j < PARAM_NUM + 1; j++) {
            for (int k = 0; k < PARAM_NUM + 1; k++) {
                eval->quadratic_mat[i][j][k] =0;
            }
        }
    }
    
    for (idrc_t i = NONLINEAR_NUM; i < zddp->size_of_polys; i++) {

        zdd_to_column(eval, i-NONLINEAR_NUM, zddp->zdd, zddp->polys[i], zddp->old_to_new);
    }
    
//    base_mat_print(eval);
    
//    printf("size_of_polys:%d\n",zddp->size_of_polys);
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
    
    
    unsigned long base[3] = {0,0,0};
    unsigned long constant[3] = {0,0,0};


    
    
    int row = 0;
    for (int i = 0; i< PARAM_NUM + VAR_NUM + 1; i++) {
        base[0] = 0;
        base[1] = 0;
        base[2] = 0;
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
                        int len = PARAM_NUM + VAR_NUM  - old_to_new->values[t->index];
                        int bn = len / 64;
                        int pos = len - bn * 64;
//                        base ^= ((unsigned long)1 << (PARAM_NUM + VAR_NUM - 1 - old_to_new->values[t->index]));
                        base[bn] ^= ((unsigned long)1 << pos);
                    }else{  //Ci's
                        
                        eval->quadratic_mat[col][row][old_to_new->values[t->index]] =1;
                    }
                    
                    e = cuddE(t);
                    while (e != zero && e != one) {
                                                      
                        if(e->index > 100){
                            int len = PARAM_NUM + VAR_NUM  -old_to_new->values[e->index];
                            int bn = len / 64;
                            int pos = len - bn * 64;
//                            base ^= ((unsigned long)1<< (PARAM_NUM + VAR_NUM - 1 -old_to_new->values[e->index]));
                            base[bn] ^= ((unsigned long)1 << pos);

                            
                            
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
                
                
                eval->base_mat[row][col][0] = base[0];
                eval->base_mat[row][col][1] = base[1];
                eval->base_mat[row][col][2] = base[2];

                
                
                
                
            }else{  //Xi*(Xi's)
                
                if(cuddT(p)!= one){
                    printf("%d not a linear polys after evaluation!!!!\n\n", col);
                    
                    break;
                }else{
                    int len = PARAM_NUM + VAR_NUM  -old_to_new->values[p->index];
                    int cn = len / 64;
                    int pos = len - cn * 64;
                    constant[cn] ^= ((unsigned long) 1 << pos);
//                    constant ^= ((unsigned long)1<< (PARAM_NUM + VAR_NUM - 1 -old_to_new->values[p->index]));
                }
                
                
                
                
                
            }
            p = cuddE(p);
            
        }
        eval->constant_array[col][0] = constant[0];
        eval->constant_array[col][1] = constant[1];
        eval->constant_array[col][2] = constant[2];

    }
    
}

void base_mat_print(evaluation_t *eval){
    
    printf("\n------------ base  matrix ---------------\n");
    
    for (int i = 0; i< PARAM_NUM; i++) {
        for (int j = 0; j < LINEAR_NUM; j++) {
            
            printf("%lu ",eval->base_mat[i][j][0]);
        }
        printf("\n");
    }
    
    printf("\n------------ constant array ---------------\n");
    
    for (int i = 0; i < LINEAR_NUM; i++) {
        printf("%lu\n",eval->constant_array[i][0] );
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
    
//    printf("\n------------ nonlinear matrices ---------------\n");
//    for (int i = 0; i < 80 - LINEAR_NUM; i++) {
//        printf("%d nonlinear polys\n",i);
//        for (int j = 0; j < PARAM_NUM+ VAR_NUM+1; j++) {
//            for (int k = 0; k < PARAM_NUM +VAR_NUM +1; k++) {
//                printf("%d ",eval->nonlinear_mat[i][j][k]);
//            }
//            printf("\n");
//        }
//    }

    
}

/*
void firsteval(evaluation_t *eval, long value, temp_t *temp){
    //B*v+c
    for (int i = 0; i < LINEAR_NUM; i++) {
        unsigned int v = 0;
        for (int j = 0; j < PARAM_NUM; j++) {
            int c = value & (1<<j);
            v ^= c*(eval->base_mat[i][PARAM_NUM - j]);
            
        }
        temp->coff_mat[i] = v^eval->constant_array[i];
        
    }
    
    //e conpute Ci*Cjs...
    
    for (int i = 0; i < LINEAR_NUM; i++) {
        for (int j = 0; j < PARAM_NUM; j++) {
            
            bool ecols[PARAM_NUM+1];
            for (int k = 0; k < PARAM_NUM; k++) {
                bool c = value & (1<<k);
                ecols[j] ^= c*eval->quadratic_mat[i][j][k];
            }
            
            ecols[j] ^= eval->quadratic_mat[i][j][PARAM_NUM];
            
            temp->constant_array[i] ^= (value & (1<<j))*ecols[j];
            
        }
        
        temp->constant_array[i] ^= eval->quadratic_mat[i][PARAM_NUM][PARAM_NUM];
    }
    
    
    
    
    
    
    
}

*/
