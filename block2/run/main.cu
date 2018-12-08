#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

//#include <helper_cuda.h>

#define DEBUG 0


#define ENUM_NUM 16        // the number of loops in each thread
#define UNKNOWN_NUM 187    // the number of unknowns
#define POLY_NUM 190       // the number of linear polynomials
#define PARA_NUM 37        // the number of parameters
#define CHECK_NUM 9       
#define SOL_MAX_NUM 200


//for GPU
#define BLOCK_NUM 32 //2^5
#define THREAD_NUM  256 // 2^8
#define THREADS_SHIFT 13 // (5+8)


//for keccak
#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define KeccakP1600_stateSizeInBytes    200

typedef long value_t; // to save values of variables.
typedef unsigned long linearpart_t; // to save 32 unknowns and 1 contants.
typedef unsigned long squarepart_t;
typedef unsigned long oripoly_t;
typedef unsigned char UINT8;
typedef unsigned long long UINT64;
typedef UINT64 tKeccakLane;


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

//Operations of keccak. On CPU

static tKeccakLane KeccakRoundConstants[maxNrRounds];
static unsigned int KeccakRhoOffsets[nrLanes];

__constant__ tKeccakLane const_KeccakRoundConstants[maxNrRounds] =
{
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
};


__constant__ unsigned int const_KeccakRhoOffsets[nrLanes] =
{
     0,  1, 62, 28, 27, 36, 44,  6, 55, 20,  3, 10, 43, 25, 39, 41, 45, 15, 21,  8, 18,  2, 61, 56, 14
};


void KeccakP1600_InitializeRoundConstants(void);
void KeccakP1600_InitializeRhoOffsets(void);
static int LFSR86540(UINT8 *LFSR);
int displayLevel = 10;


 void KeccakP1600_StaticInitialize(void)
{
    if (sizeof(tKeccakLane) != 8) {
        printf("tKeccakLane should be 64-bit wide\n");
        abort();
    }
    KeccakP1600_InitializeRoundConstants();
    KeccakP1600_InitializeRhoOffsets();
}


 void KeccakP1600_InitializeRoundConstants(void)
{
    UINT8 LFSRstate = 0x01;
    unsigned int i, j, bitPosition;

    for(i=0; i<maxNrRounds; i++) {
        KeccakRoundConstants[i] = 0;
        for(j=0; j<7; j++) {
            bitPosition = (1<<j)-1; /* 2^j-1 */
            if (LFSR86540(&LFSRstate))
                KeccakRoundConstants[i] ^= (tKeccakLane)1<<bitPosition;
        }
    }
}


 void KeccakP1600_InitializeRhoOffsets(void)
{
    unsigned int x, y, t, newX, newY;

    KeccakRhoOffsets[index(0, 0)] = 0;
    x = 1;
    y = 0;
    for(t=0; t<24; t++) {
        KeccakRhoOffsets[index(x, y)] = ((t+1)*(t+2)/2) % 64;
        newX = (0*x+1*y) % 5;
        newY = (2*x+3*y) % 5;
        x = newX;
        y = newY;
    }
}


static int LFSR86540(UINT8 *LFSR)
{
    int result = ((*LFSR) & 0x01) != 0;
    if (((*LFSR) & 0x80) != 0)
    /* Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1 */
        (*LFSR) = ((*LFSR) << 1) ^ 0x71;
    else
        (*LFSR) <<= 1;
    return result;
}

//Operations of keccak. on CPU and GPU
__host__ __device__ void KeccakP1600_Initialize(void *state)
{
    memset(state, 0, 1600/8);
}


/* ---------------------------------------------------------------- */

__host__ __device__ void KeccakP1600_AddByte(void *state, unsigned char byte, unsigned int offset)
{
    assert(offset < 200);
    ((unsigned char *)state)[offset] ^= byte;
}

/* ---------------------------------------------------------------- */

__host__ __device__ void KeccakP1600_AddBytes(void *state, const unsigned char *data, unsigned int offset, unsigned int length)
{
    unsigned int i;

    assert(offset < 200);
    assert(offset+length <= 200);
    for(i=0; i<length; i++)
        ((unsigned char *)state)[offset+i] ^= data[i];
}

/* ---------------------------------------------------------------- */

__host__ __device__ void KeccakP1600_OverwriteBytes(void *state, const unsigned char *data, unsigned int offset, unsigned int length)
{
    assert(offset < 200);
    assert(offset+length <= 200);
    memcpy((unsigned char*)state+offset, data, length);
}

/* ---------------------------------------------------------------- */

__host__ __device__ void KeccakP1600_OverwriteWithZeroes(void *state, unsigned int byteCount)
{
    assert(byteCount <= 200);
    memset(state, 0, byteCount);
}


#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)

static __host__ __device__ void theta(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane C[5]={0,0,0,0,0}, D[5]={0,0,0,0,0};

    for(x=0; x<5; x++) {
        C[x] = 0;
        for(y=0; y<5; y++)
            C[x] ^= A[index(x, y)];
    }
    for(x=0; x<5; x++)
        D[x] = ROL64(C[(x+1)%5], 1) ^ C[(x+4)%5];
    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(x, y)] ^= D[x];
}


static void rho(tKeccakLane *A)
{
    unsigned int x, y;

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(x, y)] = ROL64(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

static __device__ void rho_Device(tKeccakLane *A)
{
    unsigned int x, y;

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(x, y)] = ROL64(A[index(x, y)], const_KeccakRhoOffsets[index(x, y)]);
}

static __host__ __device__ void pi(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane tempA[25];

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        tempA[index(x, y)] = A[index(x, y)];
    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(0*x+1*y, 2*x+3*y)] = tempA[index(x, y)];
}

static __host__ __device__ void chi(tKeccakLane *A)
{
    unsigned int x, y;
    tKeccakLane C[5];

    for(y=0; y<5; y++) {
        for(x=0; x<5; x++)
            C[x] = A[index(x, y)] ^ ((~A[index(x+1, y)]) & A[index(x+2, y)]);
        for(x=0; x<5; x++)
            A[index(x, y)] = C[x];
    }
}

static void iota(tKeccakLane *A, unsigned int indexRound)
{
    A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}

static __device__ void iota_Device(tKeccakLane *A, unsigned int indexRound)
{
    A[index(0, 0)] ^= const_KeccakRoundConstants[indexRound];
}

void KeccakP1600Round(tKeccakLane *state, unsigned int indexRound)
{
#ifdef KeccakReference
    printf("11\n");
    displayRoundNumber(3, indexRound);
#endif

    theta(state);

#ifdef KeccakReference
    displayStateAsLanes(3, "After theta", state, 1600);
#endif

    rho(state);
#ifdef KeccakReference
    displayStateAsLanes(3, "After rho", state, 1600);
#endif

    pi(state);
#ifdef KeccakReference
    displayStateAsLanes(3, "After pi", state, 1600);
#endif

    chi(state);
#ifdef KeccakReference
    displayStateAsLanes(3, "After chi", state, 1600);
#endif

    iota(state, indexRound);
#ifdef KeccakReference
    displayStateAsLanes(3, "After iota", state, 1600);
#endif
}



__device__ void KeccakP1600Round_Device(tKeccakLane *state,unsigned int indexRound) {

	theta(state);
	rho_Device(state);
	pi(state);
	chi(state);
	iota_Device(state, indexRound);

}


/**
 * Initial the state from the output of block 1
 */
void stateInit(tKeccakLane state[25]) {

	KeccakP1600_StaticInitialize();
	for(int i = 0; i < 5; i++){
		state[i] = 0;
		state[i+10]=0;
	}
        state[18] = 0x90568729AB3FF556;
	state[19] = 0x16D513824AA92784;
	state[20] = 0x04A0F3940480011D;
	state[21] = 0x0A87CA1C3B68EF23;
	state[22] = 0xF5A858BD84F680F6;
	state[23] = 0x63A9789654C00AA9;
	state[24] = 0xE128EC7DB556D87B;
	for (int i = 0; i < 5; i++) {
		state[5 + i] = ~(state[20 + i]);
		if (i < 3) {
			state[15 + i] = ~(state[20 + i]);
		}
	}


}


/*
 * Obtain the initial state from the solutions 
 */

void getStates(tKeccakLane state[25], oripoly_t var_all[640][4], value_t val,
		value_t solutions[3]) {

	value_t val_sol[4];
	val_sol[3] = solutions[2];
	val_sol[2] = solutions[1];
	val_sol[1] = solutions[0];
	val_sol[0] = val ^ ((value_t) 1 << 37);

	for (int i = 0; i < 640; i++) {
		value_t w[4] = { 0, 0, 0, 0 };
		for (int j = 0; j < 4; j++) {
			w[j] = var_all[i][j] & val_sol[j];

		}

		w[0] = w[0] ^ w[1] ^ w[2] ^ w[3];
		w[0] = (w[0]) ^ (w[0] >> 32);
		w[0] = (w[0]) ^ (w[0] >> 16);
		w[0] = (w[0]) ^ (w[0] >> 8);
		w[0] = (w[0]) ^ (w[0] >> 4);
		w[0] = (w[0]) ^ (w[0] >> 2);
		w[0] = (w[0]) ^ (w[0] >> 1);
		if (w[0] & (value_t) 1) {
			int n = (i / 64 > 4) ? (i / 64 + 5) : i / 64;
			state[n] ^= ((UINT64) 1) << (i % 64);

		}
	}


}

/**
 *  Verify solutions
 */

int checkHashValue(tKeccakLane state[25]) {


	tKeccakLane state_copy[25];
	for(int i = 0; i < 25; i++){
		state_copy[i] = state[i];
	}





	for (int i = 0; i < 3; i++) {
		KeccakP1600Round(state, i);
	}



	int result = 0;

    	if (state[0] == 0xF4FE7CCEA5D8B144 && state[1] == 0x60F6C316572983A8 && state[2] == 0xA2564CA289E5F897 && ((state[3] ^= 0x00000000CA30DB85) & (0x00000000FFFFFFFF)) == 0) {

		printf("Find Preimage!!!\nState after XORed with block2:");
		for (int i = 0; i < 25; i++) {
			binary_print(state_copy[i], 64);
			printf(" ");
			printf("%llu ",state_copy[i]);
		}
		printf("\n");


		result = 1;
	}
	return result;
}


//the output of block 1 used on GPU
__constant__ tKeccakLane const_state[25] = {0x0L,0x0L,0x0L,0x0L,0x0L,
        0xFB5F0C6BFB7FFEE2L,0xF57835E3C49710DCL,0x0A57A7427B097F09L,0x9C568769AB3FF556L,0x1ED713824AA92784L,
        0x0L,0x0L,0x0L,0x0L,0x0L,
        0xFB5F0C6BFB7FFEE2L,0xF57835E3C49710DCL,0x0A57A7427B097F09L,0x90568729AB3FF556L,0x16D513824AA92784L,
        0x04A0F3940480011DL,0x0A87CA1C3B68EF23L,0xF5A858BD84F680F6L,0x63A9789654C00AA9L,0xE128EC7DB556D87BL};


__device__ linearpart_t d_linear_mat[ENUM_NUM * POLY_NUM * 3];
__device__ squarepart_t d_square_mat[ENUM_NUM * POLY_NUM];
__device__ value_t d_var_all[2560];


/**
 * Find the position of the first nonzero bit, both on GPU and CPU 
 */
static inline __host__ __device__ int largestpos(value_t val, int len) {

	for (int i = len - 1; i >= 0; i--) {
		if (val & ((value_t) 1 << i)) {
			return i;
		}
	}

	return -1;
}

/**
 * Find the position of the first nonzero bit, both on GPU and CPU 
 */
static inline __host__ __device__ int largestpos_3(value_t val0, value_t val1,
		value_t val2, int len) {
	int p = 0;
	if (len > 128) {
		p = largestpos(val2, len - 128);
		if (p > -1) {
			return p + 128;
		} else {
			p = largestpos(val1, 64);
			if (p > -1) {
				return p + 64;
			} else {
				p = largestpos(val0, 64);

			}
		}
	} else if (len > 64 && len <= 128) {
		p = largestpos(val1, len - 64);
		if (p > -1) {
			return p + 64;
		} else {
			p = largestpos(val0, 64);

		}
	} else {
		p = largestpos(val0, 64);

	}

	return p;
}

/**
 * Solve the linear system by Guassian Elimination, on CPU
 */
static inline value_t gauss_host(linearpart_t working_mat[POLY_NUM][3],
		const int poly_num, const int unknown_num, value_t solutions[SOL_MAX_NUM][3]) {

	int pos_arr[POLY_NUM]; 
	int rank = 0;

	for (int pi = 0; pi < POLY_NUM; pi++) {

		if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0 && working_mat[pi][2] == 0) {
			continue;
		}


		pos_arr[pi] = largestpos_3(working_mat[pi][0],working_mat[pi][1],working_mat[pi][2], unknown_num + 1);

		rank++;
		if (pos_arr[pi] == 0) {
			return 0;
		}





		for (int j = pi + 1; j < POLY_NUM; j++) {

			if(working_mat[j][pos_arr[pi]/64] & ((linearpart_t)1 << (pos_arr[pi] % 64))){
                working_mat[j][0] ^= (working_mat[pi][0]);
                working_mat[j][1] ^= (working_mat[pi][1]);
                working_mat[j][2] ^= (working_mat[pi][2]);
            }
		}



	}


	// back reduced
	for (int pi = 0; pi < POLY_NUM; pi++) {

		if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
				&& working_mat[pi][2] == 0) {
			continue;
		}

		for (int j = 0; j < pi; j++) {
			if (working_mat[j][pos_arr[pi] / 64]
					& ((linearpart_t) 1 << (pos_arr[pi] % 64))) {
				working_mat[j][0] ^= (working_mat[pi][0]);
				working_mat[j][1] ^= (working_mat[pi][1]);
				working_mat[j][2] ^= (working_mat[pi][2]);
			}
		}
	}

	if (rank == unknown_num) {

		// only one solution.
		solutions[0][0] = 0;
		solutions[0][1] = 0;
		solutions[0][2] = 0;
		for (int pi = 0; pi < POLY_NUM; pi++) {

			if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0 &&  working_mat[pi][2] == 0) {
				continue;
			}
			if (working_mat[pi][0] & (linearpart_t)1) {
			    solutions[0][(pos_arr[pi]-1) /64 ] ^= ((value_t)1 << (pos_arr[pi]-1) % 64);
			}
		}

		return 1;

	} else {

		// multi-solutions
		solutions[0][0] = 0;
		solutions[0][1] = 0;
		solutions[0][2] = 0;
		value_t sol_num = 1;
		bool appear[UNKNOWN_NUM + 1] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};



		for (int pi = 0; pi < POLY_NUM; pi++) {

			if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0 && working_mat[pi][2] == 0) {
				continue;
			}

			appear[pos_arr[pi]] = true;
			if (working_mat[pi][0] & (linearpart_t)1) {
			    solutions[0][(pos_arr[pi]-1) /64 ] ^= ((value_t)1 << (pos_arr[pi]-1) % 64);
			}
		}

		// duplicate solutions.
		for (int i = 1; i < UNKNOWN_NUM+1; i++) {  

			if (appear[i] == false) {


				for (int j = 0; j < sol_num; j++) {
				    
					solutions[j + sol_num][0] = solutions[j][0];
					solutions[j + sol_num][1] = solutions[j][1];
					solutions[j + sol_num][2] = solutions[j][2];
					solutions[j + sol_num][(i-1)/64] ^= ((value_t)1 << ((i-1)%64));
				}


                for (int pi = 0; pi < POLY_NUM; pi++) {
				    if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0 && working_mat[pi][2] == 0) {
			            continue;
		            }

				    for (int j = 0; j < sol_num * ((working_mat[pi][i/64] & (((linearpart_t) 1) << (i%64))) != 0); j++) {

				    	solutions[j + sol_num][(pos_arr[pi] - 1)/64] ^= ((value_t) 1 << ((pos_arr[pi] - 1)% 64));
				    }


				}


				sol_num *= 2;

			}
		}

		return sol_num;

	}

}



/**
 * Solve the linear system by Guassian Elimination, on GPU
 */
static inline __device__ value_t gauss(value_t solutions[SOL_MAX_NUM][3], linearpart_t working_mat[POLY_NUM][3],
		const int poly_num, const int unknown_num) {

	// bear revised
	int pos_arr[POLY_NUM]; // bear revised
	int rank = 0;

	for (int pi = 0; pi < POLY_NUM; pi++) {

		if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
				&& working_mat[pi][2] == 0) {
			continue;
		}

		pos_arr[pi] = largestpos_3(working_mat[pi][0], working_mat[pi][1],working_mat[pi][2],unknown_num + 1);
		rank++;


		if (pos_arr[pi] == 0) {
			return 0;
		}

		for (int j = pi + 1; j < POLY_NUM; j++) {

			if (working_mat[j][pos_arr[pi] / 64]
					& ((linearpart_t) 1 << (pos_arr[pi] % 64))) {
				working_mat[j][0] ^= (working_mat[pi][0]);
				working_mat[j][1] ^= (working_mat[pi][1]);
				working_mat[j][2] ^= (working_mat[pi][2]);
			}
		}

	}

	// back
	for (int pi = 0; pi < POLY_NUM; pi++) {

		if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
				&& working_mat[pi][2] == 0) {
			continue;
		}

		for (int j = 0; j < pi; j++) {
			if (working_mat[j][pos_arr[pi] / 64]
					& ((linearpart_t) 1 << (pos_arr[pi] % 64))) {
				working_mat[j][0] ^= (working_mat[pi][0]);
				working_mat[j][1] ^= (working_mat[pi][1]);
				working_mat[j][2] ^= (working_mat[pi][2]);
			}
		}
	}

	if (rank == unknown_num) {

		// only one solution.
		solutions[0][0] = 0;
		solutions[0][1] = 0;
		solutions[0][2] = 0;
		for (int pi = 0; pi < POLY_NUM; pi++) {

			if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
					&& working_mat[pi][2] == 0) {
				continue;
			}
			if (working_mat[pi][0] & (linearpart_t) 1) {
				solutions[0][(pos_arr[pi] - 1) / 64] ^= ((value_t) 1
						<< (pos_arr[pi] - 1) % 64);
			}
		}

		return 1;

	} else {


		solutions[0][0] = 0;
		solutions[0][1] = 0;
		solutions[0][2] = 0;
		value_t sol_num = 1;
	
		bool appear[UNKNOWN_NUM + 1] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0 };

		for (int pi = 0; pi < POLY_NUM; pi++) {

			if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
					&& working_mat[pi][2] == 0) {
				continue;
			}

			appear[pos_arr[pi]] = true;
			if (working_mat[pi][0] & (linearpart_t) 1) {
				solutions[0][(pos_arr[pi] - 1) / 64] ^= ((value_t) 1
						<< (pos_arr[pi] - 1) % 64);
			}
		}

		// duplicate solutions.
		for (int i = 1; i < UNKNOWN_NUM + 1; i++) {  
			if (appear[i] == false) {

				for (int j = 0; j < sol_num; j++) {

					solutions[j + sol_num][0] = solutions[j][0];
					solutions[j + sol_num][1] = solutions[j][1];
					solutions[j + sol_num][2] = solutions[j][2];
					solutions[j + sol_num][(i - 1) / 64] ^= ((value_t) 1
							<< ((i - 1) % 64));
				}


				for (int pi = 0; pi < POLY_NUM; pi++) {
					if (working_mat[pi][0] == 0 && working_mat[pi][1] == 0
							&& working_mat[pi][2] == 0) {
						continue;
					}

					for (int j = 0;j< sol_num* ((working_mat[pi][i / 64]& (((linearpart_t) 1)<< (i % 64))) != 0);j++) {

						solutions[j + sol_num][(pos_arr[pi] - 1) / 64] ^=((value_t) 1 << ((pos_arr[pi] - 1) % 64));
					}

				}

				sol_num *= 2;

			}
		}

		return sol_num;

	}

}




/**
 * Traverse the values of parameters to obtain matrices of linear system. Solve the system and verify the solutions.
 */
__global__ void solveLinear(const linearpart_t *d_working_mat_copy,
		const squarepart_t *d_const_mat, value_t *d_val, value_t *d_sol_total,value_t* result) {

	int thidx = blockDim.x * blockIdx.x + threadIdx.x;
	value_t val = d_val[thidx];

	linearpart_t working_mat[POLY_NUM][3]; // initialized as the const part of linear matrix. also used as the results of linear part.
	linearpart_t working_mat_copy[POLY_NUM][3];

	squarepart_t const_mat[POLY_NUM];
	d_sol_total[thidx] = 0;



	//copy data from device
	for (int i = 0; i < POLY_NUM; i++) {
		working_mat_copy[i][0] = d_working_mat_copy[thidx * POLY_NUM * 3 + i*3];
		working_mat_copy[i][1] = d_working_mat_copy[thidx * POLY_NUM * 3 + i*3 + 1];
		working_mat_copy[i][2] = d_working_mat_copy[thidx * POLY_NUM * 3 + i*3 + 2];

		const_mat[i] = d_const_mat[thidx * POLY_NUM + i];
	}

	// main loop.
	for (value_t count = 1; count < (1 << ENUM_NUM); count++) {

		// generate the next gray code
		int pos = 64-__ffsll(__brevll(count ^ (count - 1)));

		val = val ^ ((value_t) 1 << pos);


		for (int pi = 0; pi < POLY_NUM; pi++) {
			working_mat_copy[pi][0] ^= d_linear_mat[pos * POLY_NUM * 3 + pi * 3];
			working_mat_copy[pi][1] ^= d_linear_mat[pos * POLY_NUM * 3 + pi * 3 + 1];
			working_mat_copy[pi][2] ^= d_linear_mat[pos * POLY_NUM * 3 + pi * 3 + 2];

			const_mat[pi] ^= d_square_mat[pos * POLY_NUM + pi];

			working_mat[pi][0] = working_mat_copy[pi][0];
			working_mat[pi][1] = working_mat_copy[pi][1];
			working_mat[pi][2] = working_mat_copy[pi][2];


			value_t w = const_mat[pi] & val;


			working_mat[pi][0] ^= (bool)((__popcll((unsigned long long int)w)) & (value_t) 1);


		}

		value_t solutions[SOL_MAX_NUM][3];
		value_t sol_num = 0;


		// gaussian elimination
		sol_num = gauss(solutions, working_mat, POLY_NUM, UNKNOWN_NUM);
		d_sol_total[thidx] += sol_num;

		// verify on 3 round keccak.
		tKeccakLane dState[25];


		for(int s = 0;s < sol_num;s++){
			dState[0] = 0;
			dState[1] = 0;
			dState[2] = 0;
			dState[3] = 0;
			dState[4] = 0;
			dState[5] = const_state[5];
			dState[6] = const_state[6];
			dState[7] = const_state[7];
			dState[8] = const_state[8];
			dState[9] = const_state[9];
			dState[10] = 0;
			dState[11] = 0;
			dState[12] = 0;
			dState[13] = 0;
			dState[14] = 0;
			dState[15] = const_state[15];
			dState[16] = const_state[16];
			dState[17] = const_state[17];
			dState[18] = const_state[18];
			dState[19] = const_state[19];

			dState[20] = const_state[20];
			dState[21] = const_state[21];
			dState[22] = const_state[22];
			dState[23] = const_state[23];
			dState[24] = const_state[24];


			value_t val_sol[4];
			val_sol[3] = solutions[s][2];
			val_sol[2] = solutions[s][1];
			val_sol[1] = solutions[s][0];
			val_sol[0] = val ^ ((value_t)1 << 37);

			for(int i = 0; i< 640; i ++){
				value_t w[4] ={0,0,0,0};
				for(int j = 0; j< 4; j++){
					w[j] = d_var_all[i * 4 + j] & val_sol[j];
				}

				w[0] = w[0] ^w[1]^w[2]^w[3];


				if ((bool)((__popcll((unsigned long long int)w[0])) & (value_t) 1)) {
					int n = (i/64 > 4 )?( i/64 + 5 ): i/64 ;
					dState[n] ^= ((UINT64)1) << (i % 64);
				}
			}


			tKeccakLane state_copy[25];
			for(int i = 0; i < 25; i++){
				state_copy[i] = dState[i];
			}



            KeccakP1600Round_Device(dState, 0);

            KeccakP1600Round_Device(dState, 1);
            KeccakP1600Round_Device(dState, 2);



            if (dState[0] == 0xF4FE7CCEA5D8B144 && dState[1] == 0x60F6C316572983A8 && dState[2] == 0xA2564CA289E5F897 && ((dState[3] ^= 0x00000000CA30DB85) & (0x00000000FFFFFFFF)) == 0) {
				printf("Find Preimage!!! val is %lu.\n",val);
				result[0] = val;
				result[1] = val_sol[1];
				result[2] = val_sol[2];
				result[3] = val_sol[3];

				for (int i = 0; i < 25; i++) {
					printf("%llx ", state_copy[i]);

				}
				printf("\n");
			}

		}
	}

}






int main(int argc, char** argv) {
	tKeccakLane state[25];
	stateInit(state);

	const int para_num = PARA_NUM;
	const int enum_num = ENUM_NUM;


	value_t set_val = atol(argv[1])<<THREADS_SHIFT;

	const int unknown_num = UNKNOWN_NUM;
	const int poly_num = POLY_NUM;

	linearpart_t linear_mat[para_num][poly_num][3];
	linearpart_t working_mat[poly_num][3]; // initialized as the const part of linear matrix. also used as the results of linear part.
	linearpart_t working_mat_copy[poly_num][3];
	linearpart_t working_mat_file[poly_num][3];

	squarepart_t square_mat[para_num][poly_num];
	squarepart_t const_mat[poly_num]; // used to compute the const part from square polys.
	oripoly_t var_all[640][4];


	cudaSetDevice(atoi(argv[2]));

	// read the matrix files
	FILE *in1 = fopen("./59/59linear.txt", "r+");
	FILE *in2 = fopen("./59/59square.txt", "r+");
	FILE *in4 = fopen("./59/59working.txt", "r+");
	FILE *in5 = fopen("./59/59totalLinear640.txt", "r+");

	char c1, c2, c4, c5;
	for (int i = 0; i < para_num; i++) {
		for (int j = 0; j < poly_num; j++) {
			linear_mat[i][j][0] = 0;
			linear_mat[i][j][1] = 0;
			linear_mat[i][j][2] = 0;
			square_mat[i][j] = 0;

			for (int k = 0; k < 192; k++) {
				fscanf(in1, "%c", &c1);
				while (c1 != '0' && c1 != '1') {
					fscanf(in1, "%c", &c1);
				}
				if (c1 == '1') {

					linear_mat[i][j][k / 64] ^= ((linearpart_t) 1 << (k % 64));
				}
			}

			for (int k = 0; k < 64; k++) {
				fscanf(in2, "%c", &c2);
				while (c2 != '0' && c2 != '1') {
					fscanf(in2, "%c", &c2);
				}
				if (c2 == '1') {
					square_mat[i][j] ^=
							((squarepart_t) 1 << (para_num - 1 - k));
				}
			}
		}

	}

	for (int i = 0; i < poly_num; i++) {
		working_mat[i][0] = 0;
		working_mat[i][1] = 0;
		working_mat[i][2] = 0;
		for (int j = 0; j < 192; j++) {
			fscanf(in4, "%c", &c4);
			while (c4 != '0' && c4 != '1') {
				fscanf(in4, "%c", &c4);
			}
			if (c4 == '1') {

				working_mat[i][j / 64] ^= ((linearpart_t) 1 << (j % 64));
			}
		}
		working_mat_file[i][0] = working_mat[i][0];
		working_mat_file[i][1] = working_mat[i][1];
		working_mat_file[i][2] = working_mat[i][2];
	}

	for (int i = 0; i < 640; i++) {
		var_all[i][0] = 0;
		var_all[i][1] = 0;
		var_all[i][2] = 0;
		var_all[i][3] = 0;
		for (int j = 0; j < 256; j++) {

			fscanf(in5, "%c", &c5);
			while (c5 != '0' && c5 != '1') {
				fscanf(in5, "%c", &c5);

			}
			if (c5 == '1') {
				var_all[i][j / 64] ^= ((value_t) 1 << (j % 64));
			}

		}

	}

	fclose(in1);
	fclose(in2);
	fclose(in4);
	fclose(in5);

	printf("finish reading file!\n");


	//allocate device memory

	linearpart_t linear_mat_enum[ENUM_NUM * POLY_NUM * 3];
	squarepart_t square_mat_enum[ENUM_NUM * POLY_NUM];
	value_t var_all_enum[640 * 4];

	for (int i = 0; i < ENUM_NUM; i++) {
		for (int j = 0; j < POLY_NUM; j++) {
			for (int k = 0; k < 3; k++) {
				linear_mat_enum[i * POLY_NUM * 3 + j * 3 + k] =
						linear_mat[i][j][k];

			}
		}
	}

	for (int i = 0; i < ENUM_NUM; i++) {
		for (int j = 0; j < POLY_NUM; j++) {
			square_mat_enum[i * POLY_NUM + j] = square_mat[i][j];
		}
	}

	cudaMemcpyToSymbol(d_linear_mat, linear_mat_enum,
			3 * ENUM_NUM * POLY_NUM * sizeof(linearpart_t));
	cudaMemcpyToSymbol(d_square_mat, square_mat_enum,
			ENUM_NUM * POLY_NUM * sizeof(squarepart_t));

	for (int i = 0; i < 640; i++) {
		for (int j = 0; j < 4; j++) {
			var_all_enum[i * 4 + j] = var_all[i][j];
		}

	}
	cudaMemcpyToSymbol(d_var_all, var_all_enum, 640 * 4 * sizeof(value_t));

	printf("finish copying device memory!\n");

	cudaError_t err = cudaSuccess;
	int thidx = BLOCK_NUM * THREAD_NUM;

	value_t *d_val = NULL;
	err = cudaMalloc((void **) &d_val, thidx * sizeof(value_t));
	if (err != cudaSuccess) {
		printf("Failed to allocate device value (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	value_t *d_sol_total = NULL;
	err = cudaMalloc((void **) &d_sol_total, thidx * 3 * sizeof(value_t));
	if (err != cudaSuccess) {
		printf("Failed to allocate device value (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	linearpart_t *d_working_mat_copy = NULL;
	err = cudaMalloc((void **) &d_working_mat_copy,
			thidx * poly_num * 3 * sizeof(linearpart_t));
	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to allocate device working_mat_copy (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	squarepart_t *d_const_mat = NULL;
	err = cudaMalloc((void **) &d_const_mat,
			thidx * poly_num * sizeof(squarepart_t));
	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to allocate devices const_mat (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to copy oripolys from host to device (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	value_t h_result[4] = {0, 0, 0, 0};
	value_t *d_result = NULL;
	err = cudaMalloc((void **) &d_result, 4 *sizeof(value_t));
	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to allocate devices result (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	err = cudaMemcpy(d_result, h_result,
	4 * sizeof(value_t),cudaMemcpyHostToDevice);

		if (err != cudaSuccess) {
			fprintf(stderr,
					"Failed to copy result from host to device (error code %s)!\n",
					cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}


	printf("finish allocate device memory!\n");


	//deal with the case with all parameters are 0's on CPU
	int s_total_p0 = 0;
	value_t *val_arr = (value_t*)calloc(thidx, sizeof(value_t));
	linearpart_t *working_mat_copy_arr = (linearpart_t*)calloc(thidx * POLY_NUM * 3, sizeof(linearpart_t));
	squarepart_t *const_mat_arr = (squarepart_t*)calloc(thidx * POLY_NUM, sizeof(squarepart_t));
	for (int thi = 0; thi < thidx; thi++) {

		value_t sol_num = 0;
		value_t solutions[SOL_MAX_NUM][3];

		//int sol_total = 0;
		value_t val = (set_val + (value_t) thi) << ENUM_NUM;

		val_arr[thi] = val;
		for (int pi = 0; pi < POLY_NUM; pi++) {
			working_mat[pi][0] = working_mat_file[pi][0];
			working_mat[pi][1] = working_mat_file[pi][1];
			working_mat[pi][2] = working_mat_file[pi][2];

			const_mat[pi] = 0;
		}

		for (int pos = enum_num; pos < para_num; pos++) {

			if (val & ((value_t) 1 << pos)) {

				for (int pi = 0; pi < poly_num; pi++) {
					working_mat[pi][0] ^= linear_mat[pos][pi][0];
					working_mat[pi][1] ^= linear_mat[pos][pi][1];
					working_mat[pi][2] ^= linear_mat[pos][pi][2];
				}

				for (int pi = 0; pi < poly_num; pi++) {
					const_mat[pi] ^= square_mat[pos][pi];

				}

			}

		}


		for (int i = 0; i < poly_num; i++) {
			working_mat_copy[i][0] = working_mat[i][0];
			working_mat_copy[i][1] = working_mat[i][1];
			working_mat_copy[i][2] = working_mat[i][2];

			working_mat_copy_arr[thi * POLY_NUM * 3 + 3 * i] = working_mat_copy[i][0];
			working_mat_copy_arr[thi * POLY_NUM * 3 + 3 * i + 1] = working_mat_copy[i][1];
			working_mat_copy_arr[thi * POLY_NUM * 3 + 3 * i + 2] = working_mat_copy[i][2];

			const_mat_arr[thi * POLY_NUM + i] = const_mat[i];

		}

		for (int pi = 0; pi < poly_num; pi++) {

			value_t w = const_mat[pi] & val;

			w = (w) ^ (w >> 32);
			w = (w) ^ (w >> 16);
			w = (w) ^ (w >> 8);
			w = (w) ^ (w >> 4);
			w = (w) ^ (w >> 2);
			w = (w) ^ (w >> 1);

			if (w & (value_t) 1) {

				working_mat[pi][0] ^= (linearpart_t) 1;
			}

		}

		sol_num = gauss_host(working_mat, POLY_NUM, UNKNOWN_NUM, solutions);
		s_total_p0 += sol_num;
		 for (int s = 0; s < sol_num; s++) {

			tKeccakLane state[25];
			stateInit(state);
			getStates(state, var_all,val, solutions[s]);

			if(checkHashValue(state)){

				FILE *out = fopen("preimage_b2_59.txt","a+");
				printf("we have done on GPU!!! val:%lu, sol:%lu %lu %lu\n",val,solutions[s][0],solutions[s][1],solutions[s][2]);
				fprintf(out,"we have done on GPU!!! val:%lu, sol:%lu %lu %lu\n",val,solutions[s][0],solutions[s][1],solutions[s][2]);
				fclose(out);
			}

		    }


	}

	printf("finish cpu computing!\n");

	//begin device part: copy value from host ot devide
	err = cudaMemcpy(d_val, val_arr, thidx * sizeof(value_t),
			cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		printf("Failed to copy value from host to device (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	err = cudaMemcpy(d_working_mat_copy, working_mat_copy_arr,
			thidx * 3 * poly_num * sizeof(linearpart_t), cudaMemcpyHostToDevice);

	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to copy working_mat_copy from host to device (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	err = cudaMemcpy(d_const_mat, const_mat_arr,
			thidx * poly_num * sizeof(squarepart_t), cudaMemcpyHostToDevice);

	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to copy const_mat from host to device (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	printf("enum num : %d\nblock num : %d\nthread num : %d\n", ENUM_NUM,
			BLOCK_NUM, THREAD_NUM);


	//solve linear system on GPU
	cudaEvent_t start1;
	cudaEventCreate(&start1);
	cudaEvent_t stop1;
	cudaEventCreate(&stop1);
	cudaEventRecord(start1, NULL);

	printf("begin solve linear system!\n");
	solveLinear<<<BLOCK_NUM, THREAD_NUM>>>(d_working_mat_copy, d_const_mat,
			d_val, d_sol_total,d_result);

	cudaEventRecord(stop1, NULL);
	cudaEventSynchronize(stop1);
	float msecTotal1 = 0.0f;
	cudaEventElapsedTime(&msecTotal1, start1, stop1);
	err = cudaGetLastError();

	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to launch solveLinear kernel (error code %s)!\n",
				cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	value_t h_sol_total[thidx];
	err = cudaMemcpy(h_sol_total, d_sol_total, thidx * sizeof(value_t),
				cudaMemcpyDeviceToHost);

	if (err != cudaSuccess) {
		fprintf(stderr,
				"Failed to copy total solution numbers from device to host (error code %s)!\n",
				cudaGetErrorString(err));
			exit(EXIT_FAILURE);
	}

	err = cudaMemcpy(h_result, d_result,4 * sizeof(value_t),
					cudaMemcpyDeviceToHost);

	if (err != cudaSuccess) {
		fprintf(stderr,
					"Failed to copy result from device to host (error code %s)!\n",
					cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

	//print the results

	if(h_result[0]!=0 || h_result[1]!=0 || h_result[2]!=0 ||  h_result[3]!=0){
	    	FILE *out = fopen("preimage_b2_59.txt","a+");
			printf("we have done on GPU!!! val:%lu, sol:%lu %lu %lu\n",h_result[0],h_result[1],h_result[2],h_result[3]);
			fprintf(out,"we have done on GPU!!! val:%lu, sol:%lu %lu %lu\n",h_result[0],h_result[1],h_result[2],h_result[3]);
			fclose(out);
		}else{
			FILE *out = fopen("result.txt","a+");
			long sol_all_threads = s_total_p0;
			for(int i = 0;i < thidx;i++){
				sol_all_threads += h_sol_total[i];
			}
			printf("val : %lu~%lu ,find %lu solutions, none is correct...\n",set_val << ENUM_NUM ,(set_val << ENUM_NUM)+(THREAD_NUM * BLOCK_NUM) * (1 << ENUM_NUM) -1, sol_all_threads);
			fprintf(out, "Part %d finished -- val : %lu~%lu ,find %lu solutions, none is correct...\n",atol(argv[1]), set_val << ENUM_NUM ,(set_val << ENUM_NUM)+(THREAD_NUM * BLOCK_NUM) * (1 << ENUM_NUM) -1, sol_all_threads);
			fclose(out);

		}


	printf("time:%.3lf ms\n---------------------------------------\n", msecTotal1);


	cudaFree(val_arr);
	cudaFree(d_working_mat_copy);
	cudaFree(d_const_mat);
	cudaFree(d_val);
	cudaFree(d_sol_total);
	cudaFree(d_result);
}
