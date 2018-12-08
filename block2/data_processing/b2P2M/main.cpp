


#include "evaluate.h"


void printBinary1(linearpart_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}


void printBinary2(squarepart_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}

void printBinary3(oripoly_t n){
    for (int i = 0; i < 64; i++) {
        printf("%lu",(n>> i) &(unsigned long)1 );
    }
}


void printBinary1_file(linearpart_t n,FILE *f){
    for (int i = 0; i < 64; i++) {
        fprintf(f,"%lu",(n>> i) &(unsigned long)1 );
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

/**
 * read meaningful char from file.
 */
static inline void meaningchar(FILE *f, char *p_c) { //r
    do {
        fscanf(f, "%c", p_c);
    } while ((*p_c) != 'x' && !((*p_c) >= '0' && (*p_c) <= '9') && (*p_c) != '(' && (*p_c) != ')' && (*p_c) != '.' && (*p_c) != ';' && (*p_c) != ',' && (*p_c) != '+');
}


/**
 * write the matrix of totalLinear file
 */
void printTotalLinearMat_file(FILE *f,FILE *out){
    char c;
    meaningchar(f, &c);
    bool TLmat[737][738];
    for (int i = 0; i < 737; i++) {
        for (int j = 0; j < 738; j++) {
            TLmat[i][j] = 0;
        }
    }
    
    bool first = true;
    int row = 0;
    while (c != '.') {
        if(c == 'x'){
            
            int d = 0;
            meaningchar(f, &c);
            while (c <= '9' && c >= '0') {
                d = d * 10 + c - '0';
                meaningchar(f, &c);
            }
            
            if(first){
                first = false;
                if(row < d){
                    row++;
                    while (row < d) {
                        TLmat[row][row] = 1;
                        row++;
                    }
                }
            }else{
                TLmat[row][d] = 1;
            }
            
        }else if (c == ','){
            first = true;
            meaningchar(f, &c);
        }else if(c == '1'){
            TLmat[row][737] = 1;
            meaningchar(f, &c);
        }else{
            meaningchar(f, &c);
        }
        
    }
    
    for (int r = row + 1; r < 737; r++) {
        TLmat[r][r] = 1;
    }
    
    
    
    int p[256] = {700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,639,638,637,636,635,634,633,632,631,630,629,628,627,626,625,624,623,622,621,620,619,618,617,616,615,614,613,612,611,610,609,608,607,606,605,604,603,602,601,600,599,598,597,596,595,594,593,592,591,590,589,588,587,586,585,584,583,582,581,580,579,578,577,575,574,573,572,571,570,569,568,567,566,565,564,563,562,561,560,559,558,557,556,555,553,552,551,550,549,548,547,546,545,543,542,541,540,539,538,537,536,535,534,533,532,531,530,529,528,527,526,525,524,523,522,521,520,519,518,517,516,515,514,513,512,511,510,509,508,507,506,505,504,503,502,501,500,499,498,497,496,495,494,493,492,491,490,489,488,487,486,483,482,481,480,479,478,477,476,475,474,473,472,471,470,469,468,467,466,465,463,462,461,460,459,458,457,456,455,454,453,452,451,450,449,448,447,0,0,0,0,0};
    
    
    for (int r = 0; r < 640; r ++) {
        for (int c =0; c < 256; c++) {
            if (c % 64 == 0) {
                fprintf(out," ");
            }
            if(p[c] > 0){
                fprintf(out,"%d", TLmat[r][p[c]]);
            }else{
                fprintf(out,"0");
            }
            
        }
        fprintf(out, "\n");
    }
    
    
}

/**
 * print the matrix of totalLinear file to screen
 */
void printTotalLinearMat(FILE *f){
    char c;
    meaningchar(f, &c);
    bool TLmat[737][738];
    for (int i = 0; i < 737; i++) {
        for (int j = 0; j < 738; j++) {
            TLmat[i][j] = 0;
        }
    }
    
    bool first = true;
    int row = 0;
    while (c != '.') {
        if(c == 'x'){
            
            int d = 0;
            meaningchar(f, &c);
            while (c <= '9' && c >= '0') {
                d = d * 10 + c - '0';
                meaningchar(f, &c);
            }
            
            if(first){
                first = false;
                if(row < d){
                    row++;
                    while (row < d) {
                        TLmat[row][row] = 1;
                        row++;
                    }
                }
            }else{
                TLmat[row][d] = 1;
            }
            
        }else if (c == ','){
            first = true;
            meaningchar(f, &c);
        }else if(c == '1'){
            TLmat[row][737] = 1;
            meaningchar(f, &c);
        }else{
            meaningchar(f, &c);
        }
        
    }
    
    for (int r = row + 1; r < 737; r++) {
        TLmat[r][r] = 1;
    }
    

    
    int p[256] = {700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,639,638,637,636,635,634,633,632,631,630,629,628,627,626,625,624,623,622,621,620,619,618,617,616,615,614,613,612,611,610,609,608,607,606,605,604,603,602,601,600,599,598,597,596,595,594,593,592,591,590,589,588,587,586,585,584,583,582,581,580,579,578,577,575,574,573,572,571,570,569,568,567,566,565,564,563,562,561,560,559,558,557,556,555,553,552,551,550,549,548,547,546,545,543,542,541,540,539,538,537,536,535,534,533,532,531,530,529,528,527,526,525,524,523,522,521,520,519,518,517,516,515,514,513,512,511,510,509,508,507,506,505,504,503,502,501,500,499,498,497,496,495,494,493,492,491,490,489,488,487,486,483,482,481,480,479,478,477,476,475,474,473,472,471,470,469,468,467,466,465,463,462,461,460,459,458,457,456,455,454,453,452,451,450,449,448,447,0,0,0,0,0};
    
    
    for (int r = 0; r < 640; r ++) {
        for (int c =0; c < 256; c++) {
            if (c % 64 == 0) {
                printf(" ");
            }
            if(p[c] > 0){
                printf("%d", TLmat[r][p[c]]);
            }else{
                printf("0");
            }
            
        }
        printf("\n");
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
    linearpart_t linear_mat[para_num][poly_num][3];
    linearpart_t working_mat[poly_num][3]; // initialized as the const part of linear matrix. also used as the results of linear part.
    linearpart_t working_mat_copy[poly_num][3];
    
    squarepart_t square_mat[para_num][poly_num];
    squarepart_t const_mat[poly_num]; // used to compute the const part from square polys.
    
    
    oripoly_t polys[ori_num][para_num + unknown_num + 1][3];
    oripoly_t cstpoly[2];
    
    
    /**
     ****************************************
     ** convert input polynomials to matrices
     ****************************************
     */
    FILE *in = fopen("../keccakEquations/59rename.txt", "r+");
    CHECK_NN(in)
    
    FILE *linear_file = fopen("../Mat/linear_mat.txt", "w+");
    FILE *working_file = fopen("../Mat/working_mat.txt", "w+");
    FILE *square_file = fopen("../Mat/square_mat.txt", "w+");
    
    
    // initialize a zdd data structure.
    zddpoly_t *zddp = zddpoly_init();
    
    zddpoly_readfile(zddp, in);
    fclose(in);
    
    
    evaluation_t eval;
    evaluation_init(&eval, zddp);

    printf("\n");
    
    // write the matrix files
    for (int i = 0; i < para_num; i++) {
    
        for (int j = 0; j < poly_num; j++) {
            linear_mat[para_num- i -1][j][0] = ((&eval)->base_mat[i][j][0])^(&eval)->quadratic_mat[j][i][para_num];
            linear_mat[para_num - i -1][j][1] = ((&eval)->base_mat[i][j][1]);
            linear_mat[para_num - i -1][j][2] = ((&eval)->base_mat[i][j][2]);
//
            printBinary1_file(linear_mat[para_num- i -1][j][0],linear_file);
            fprintf(linear_file,"|");
            printBinary1_file(linear_mat[para_num- i -1][j][1],linear_file);
            fprintf(linear_file,"|");
            printBinary1_file(linear_mat[para_num- i -1][j][2],linear_file);
            fprintf(linear_file," ");
            
        }
        fprintf(linear_file, "\n ");

    }
    
    
    
    for (int j = 0; j < poly_num; j++) {
        working_mat[j][0] = ((&eval)->constant_array[j][0])^((&eval)->quadratic_mat[j][para_num][para_num]);
        working_mat[j][1] = ((&eval)->constant_array[j][1]);
        working_mat[j][2] = ((&eval)->constant_array[j][2]);

        printBinary1_file(working_mat[j][0],working_file);
        fprintf(working_file,"|");
        printBinary1_file(working_mat[j][1],working_file);
        fprintf(working_file, "|");
        printBinary1_file(working_mat[j][2],working_file);
        fprintf(working_file,"\n");


    }
    
    
    for (int i = 0; i < para_num; i++) {
        
        for (int j = 0; j < poly_num; j++) {
            squarepart_t value = 0;
            for (int k = 0; k < para_num; k++) {
                value = (value<<1)^(squarepart_t)(&eval)->quadratic_mat[j][i][k];
            }
            square_mat[para_num-i-1][j] = value;
            printBinary2_file(square_mat[para_num-i-1][j],square_file);
            fprintf(square_file," ");
            
        }
        fprintf(square_file,"\n");
    }
    
    fclose(linear_file);
    fclose(working_file);
    fclose(square_file);
    
    FILE *fm_in = fopen("../keccakEquations/59totalLinear.txt", "r+");
    FILE *totalLinear_file = fopen("../Mat/totalLinear640.txt", "w+");
    printTotalLinearMat_file(fm_in,totalLinear_file);
    fclose(fm_in);
    fclose(totalLinear_file);
    

    exit(0);
    
}















































