



#include "evaluation.h"



void printBinary1(linearpart_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}


void printBinary1_file(linearpart_t n,FILE *f){
    for (int i = 0; i < 64; i++) {
        fprintf(f,"%lu",(n>> i) &(unsigned long)1 );
    }
}

void printBinary2(squarepart_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}

void printBinary2_file(squarepart_t n,FILE *f){
    for (int i = 0; i < 64; i++) {
        fprintf(f,"%lu",(n>>i) &(unsigned long)1 );
    }
    
}

void printBinary3_file(oripoly_t n,FILE *f){
    for (int i = 0; i < 64; i++) {
        fprintf(f,"%lu",(n>> i) &(unsigned long)1 );
    }
}

void printBinary3(oripoly_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}



int main(int argc, const char * argv[]) {
    
    
    
      // generate random matrices.
    const int para_num = PARAM_NUM;
    const int unknown_num = VAR_NUM;
    const int poly_num = LINEAR_NUM;
    const int ori_num = NONLINEAR_NUM;
    
    
    
    
    
    
    // linear part.
    
    // linear part can be unsigned int, if the unknowns are smaller than 31
    linearpart_t linear_mat[para_num][poly_num][2];
    linearpart_t working_mat[poly_num][2]; // initialized as the const part of linear matrix. also used as the results of linear part.
    linearpart_t working_mat_copy[poly_num][2];
    
    squarepart_t square_mat[para_num][poly_num];
    squarepart_t const_mat[poly_num]; // used to compute the const part from square polys.
    
    
    oripoly_t polys[ori_num][para_num + unknown_num + 1][3];
    oripoly_t cstpoly[2];
    
    
    /**
     ****************************************
     ** convert input polynomials to matrices
     ****************************************
     */
     
    FILE *linear_file = fopen("../Mat/linear_mat.txt", "w+");
    FILE *working_file = fopen("../Mat/working_mat.txt", "w+");
    FILE *square_file = fopen("../Mat/square_mat.txt", "w+");
    FILE *poly_file = fopen("../Mat/poly_mat.txt", "w+"); 
    
    FILE *in = fopen("../keccakEquations/224rename.txt", "r+");
    CHECK_NN(in)
    
    if (!in) {
     	printf("\n error in open file.\n");
     	exit(0);
    }
    
    // initialize a zdd data structure.
    zddpoly_t *zddp = zddpoly_init();
    zddpoly_readfile(zddp, in);
    
    fclose(in);
    

    evaluation_t eval;

    evaluation_init(&eval, zddp);

    printf("\n");


    //write the matrix files    
    for (int i = 0; i < para_num; i++) {
    
        for (int j = 0; j < poly_num; j++) {
            linear_mat[para_num- i -1][j][0] = ((&eval)->base_mat[i][j]<<1)^(&eval)->quadratic_mat[j][i][para_num];
            linear_mat[para_num - i -1][j][1] = ((&eval)->base_mat[i][j] >> 63);
            
            printBinary1_file(linear_mat[para_num- i -1][j][0], linear_file);
            fprintf(linear_file, "|");
            printBinary1_file(linear_mat[para_num- i -1][j][1], linear_file);
            fprintf(linear_file, " ");

            
        }

        fprintf(linear_file, "\n");

    }
    

    
    for (int j = 0; j < poly_num; j++) {
        working_mat[j][0] = ((&eval)->constant_array[j]<<1)^((&eval)->quadratic_mat[j][para_num][para_num]);
        working_mat[j][1] = ((&eval)->constant_array[j]>>63);
        printBinary1_file(working_mat[j][0],working_file);
        fprintf(working_file,"|");
        printBinary1_file(working_mat[j][1],working_file);
        fprintf(working_file,"\n");


    }
    

    
    
    for (int i = 0; i < para_num; i++) {
        
        for (int j = 0; j < poly_num; j++) {
            squarepart_t value = 0;
            for (int k = 0; k < para_num; k++) {
                value = (value<<1)^(squarepart_t)(&eval)->quadratic_mat[j][i][k];
            }
            square_mat[para_num-i-1][j] = value;
            
            printBinary2_file(square_mat[para_num-i-1][j], square_file);
            fprintf(square_file, " ");

            
        }

        fprintf(square_file, "\n");
    }
    
    for (int j = 0; j < poly_num; j++) {
        const_mat[j] = 0;
    }
    
    
    for (int i = 0; i< ori_num; i++) {
        for (int j = 0; j < para_num + unknown_num + 1; j++) {
            oripoly_t value = 0;
            // parameters
            for ( int k = 0; k < para_num ; k++) {
                value = (value<<1) ^((&eval)->nonlinear_mat[i][j][k]);

                
            }
            if(j < para_num){
                polys[i][para_num - j -1 ][0] = value;
            }else if(j == para_num + unknown_num){
                polys[i][j][0] = value;
            }else{
                polys[i][2*para_num + unknown_num -1 -j][0] = value;
            }
            
            value = 0;
            
            //unknows
            for (int k = para_num; k < para_num + unknown_num ; k++) {
                value = (value<<1) ^(&eval)->nonlinear_mat[i][j][k];

            }
            if(j < para_num){
                polys[i][para_num - j -1 ][1] = value;
            }else if(j == para_num + unknown_num){
                polys[i][j][1] = value;
            }else{
                polys[i][2*para_num + unknown_num -1 -j][1] = value;
            }

            //constants
            value = 0;
            int k = para_num + unknown_num;
            value = (&eval)->nonlinear_mat[i][j][k];
            if(j < para_num){
                polys[i][para_num - j -1 ][2] = value;
            }else if(j == para_num + unknown_num){
                polys[i][j][2] = value;
            }else{
                polys[i][2*para_num + unknown_num -1 -j][2] = value;
            }
        
        }
        
    }
    

    for (int i = 0; i< ori_num; i++) {
        for (int j = 0; j < para_num + unknown_num + 1; j++) {
            printBinary3_file(polys[i][j][0],poly_file);
            fprintf(poly_file,"|");
            printBinary3_file(polys[i][j][1],poly_file);
            fprintf(poly_file,"|");
            printBinary3_file(polys[i][j][2],poly_file);
            fprintf(poly_file," ");

        }

        fprintf(poly_file, "\n");
        
    }
    fclose(linear_file);
    fclose(working_file);
	fclose(square_file);
	fclose(poly_file);
	
	exit(0);
    
}















































