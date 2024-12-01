#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#ifdef _OPENMP 
  #include <omp.h>
#endif

#define MIN 0.0    //min value of random value
#define MAX 10.0   //max value of random value      
#define OMPTHREADS 64  //max number of omp threads 
#define NUM_RUNS 10  //number of runs to test algorithm

int n;    //matrix size
double wstart, wend;  //variable for wall clock time
double wtimeSeq, wtimeImp, wtimePar;  //safe different value of wstart - wend
double avg_time_seqCS, avg_time_seqT, avg_time_imp, avg_time_parallel, avg_speedup, avg_efficiency;  //avg values of analyses parameters 
double total_time_parallel;
double data_transferred; 
double bandwidthSeq, bandwidthImp, bandwidthPar;

double random_float(double min, double max);
void printMat (double **matrix);  //print matrix
bool checkMat (double **m1, double **m2, double **m3);  //check correctness of each matrix transpose
bool checkPowerTwo(const int n);
//sequential functions 
double** matTranspose(double **matrix);
bool checkSym(double **matrix);
//implicit parallelization functions 
double** matTransposeImp(double **matrix);
bool checkSymImp(double **matrix);
//Explicit Parallelization with OpenMP functions
double** matTransposeOMP(double **matrix);
bool checkSymOMP(double **matrix);


int main(int argc, char *argv[]){
    srand(time(NULL));

    int i, j; //indexes
    int num_threads;
    int run;
    double **M, **TS, **TI, **TP;  //matrix
    
    if (argc != 2 || !checkPowerTwo(atoi(argv[1]))) {
      	fprintf(stderr,"usage: scriptname <integer size_matrix> or <not power of 2>\n");
      	return -1;
    }
    
    //get matrix size
    n = atoi(argv[1]);
    //define amount data transfered
    data_transferred = 2.0 * n * n * sizeof(double);
    
    //allocate space for matrix
    M = (double**)malloc(n*sizeof(double*));
    
    for(i=0; i<n; i++){
        M[i] = (double*)malloc(n*sizeof(double));         
    }
    
    //Initialize random n x n matrix
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            M[i][j] = random_float(MIN, MAX);
        }
    }
    
    printf("\n\t\t________ COMPARE IMPLEMENTATIONS (matrix:%dx%d) ________\n\n",n, n);
    
    //Task 1: Sequential Implementation
    //wall clock time of checkSym
    printf("1. Sequential Implementation Times:\n");
    
    total_time_parallel = 0.0;
    for(run = 0; run < NUM_RUNS; run++){
        wstart=omp_get_wtime();
        
        checkSym(M);
        
        wend=omp_get_wtime();
        wtimeSeq = wend-wstart;
        total_time_parallel += wtimeSeq; 
    }
    
    avg_time_seqCS = total_time_parallel / NUM_RUNS; 
    
    printf("\t-checkSym function: average (on %d runs) wall clock time  = %8.4g sec\n", NUM_RUNS, avg_time_seqCS);
    
    
    //wall clock time of matTranspose
    TS = matTranspose(M);
    
    printf("\t-matTranspose function: average (on %d runs) wall clock time = %8.4g sec\n\n", NUM_RUNS, avg_time_seqT);
    
    //Task 2: Implicit Parallelization
    //wall clock time of checkSym
    printf("2. Implicit Parallelization Times:\n");
    
    total_time_parallel = 0.0;
    for(run = 0; run < NUM_RUNS; run++){
        wstart=omp_get_wtime();
        
        checkSymImp(M);
        
        wend=omp_get_wtime();
        wtimeImp = wend-wstart;
        total_time_parallel += wtimeImp; 
    }
    
    avg_time_imp = total_time_parallel / NUM_RUNS; 
    printf("\t-checkSym function: average (on %d runs) wall clock time = %8.4g sec\n", NUM_RUNS, avg_time_imp);
    
    // wall clock time of matTranspose
    TI = matTransposeImp(M);
    
    printf("\t-matTranspose function: average (on %d runs) wall clock time = %8.4g sec\n\n", NUM_RUNS, avg_time_imp);

    //Task 3: Explicit Parallelization with OpenMP
    //wall clock time of checkSym
    printf("3. Explicit Parallelization with OpenMP Average Times:\n");
    
#ifdef _OPENMP
    printf("\t-checkSymOMP function analyses (on %d runs): \n", NUM_RUNS);
    printf("Num_Threads | Avg_Parallel_Time | Avg_Speedup | Avg_Efficiency\n"); 
    // Test performance with thread counts of 1, 2, 4, 8, 16, 32, and 64
    for (int num_threads = 2; num_threads <= 64; num_threads *= 2) {
        omp_set_num_threads(num_threads);

        total_time_parallel = 0.0;

        for (int run = 0; run < NUM_RUNS; run++) {
            //start
            wstart = omp_get_wtime();

            checkSymOMP(M);

            wend = omp_get_wtime();
            wtimePar = wend - wstart;
            total_time_parallel += wtimePar;
        }

        // Compute average parallel time, speedup, and efficiency
        avg_time_parallel = total_time_parallel / NUM_RUNS;
        avg_speedup = avg_time_seqCS / avg_time_parallel;
        avg_efficiency = avg_speedup / num_threads;

        // Print results for each thread count
        printf("%11d | %17f | %11f | %9.2f%%\n", 
               num_threads, avg_time_parallel, avg_speedup, avg_efficiency * 100);
    }
    
    printf("\n\t-matTransposeOMP function analyses (on %d runs): \n", NUM_RUNS);
    // wall clock time of matTranspose
    TP = matTransposeOMP(M);

    //safe wend - wstart
    
#endif
    
    //compute effective bandwidth
    
    bandwidthSeq = data_transferred / (avg_time_seqT * 1e9); // GB/s
    bandwidthImp = data_transferred / (avg_time_imp * 1e9);  // GB/s
    
    printf("\n\t\t________ PERFORMANCE ANALYSIS ________\n\n");
    printf("3.Effective bandwdths:\n");
    printf("\t-bandwidth imp. Sequential: %lf GB/s\n", bandwidthSeq);
    printf("\t-bandwidth imp. Implicit Parallelization: %lf GB/s\n", bandwidthImp);
    
	/*
    //check if each matrix transpose goes well in each approach
    if(checkMat(TS, TI, TP)){
        printf("Matrix transpose goes well in each implemetation!\n");
    }else{
        printf("ERROR! something goes wrong some implemetation.\n");
    }*/
    
    //free allocation space of M and T
    for( i=0; i<n; i++){
        free(M[i]);
        free(TS[i]);
        free(TI[i]);
        free(TP[i]);
    }
    
    free(M);
    free(TS);
    free(TI);
    free(TP);
    
    
    return 0;
}


//
//_____________FUNCTIONS IMPLEMENTATION___________
//


void printMat (double **matrix){
    //print matrix
    int i, j;
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
          printf("%f\t", matrix[i][j]);
      }
      printf("\n");
    }
}

double random_float(double min, double max) {
    //return (float)(rand()) / (float)(rand());  //without range [min, max]
    return min + (double)rand() / RAND_MAX * (max - min);  //with range [min, max]
}

bool checkPowerTwo(const int n){
    
    if(n <= 0){
        return false;
    }
    
    return (n & (n-1)) == 0;
}

bool checkMat (double **m1, double **m2, double **m3){
    int i, j;
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(m1[i][j] != m2[i][j] && m2[i][j] != m3[i][j]){
                return false; 
            }
        }
    }
    
    return true; 
}

double** matTranspose(double **matrix){
    int i, j;
    int run;
    double **T = (double**)malloc(n*sizeof(double*));
    
    for(i=0; i<n; i++){
        T[i] = (double*)malloc(n*sizeof(double)); 
    }
    
    total_time_parallel = 0.0; 
    for(run = 0; run < NUM_RUNS; run++){
        //start time
        wstart=omp_get_wtime();
        
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                T[i][j] = matrix[j][i];
            }
        }
        //end time
        wend=omp_get_wtime();
        wtimeSeq = wend - wstart;
        total_time_parallel += wtimeSeq;
    }
    
    avg_time_seqT = total_time_parallel / NUM_RUNS;
    return T;
}

bool checkSym(double **matrix){
    int i, j;
    bool sym = true;
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(matrix[i][j] != matrix[j][i]){
                sym = false;
            }
        }
    }
    
    return sym;
}

#pragma GCC optimize("O2")
double** matTransposeImp(double **matrix){
    int i, j;
    int run;
    double **T = (double**)malloc(n*sizeof(double*));
    #pragma simd
    for(i=0; i<n; i++){
        T[i] = (double*)malloc(n*sizeof(double)); 
    }
    
    total_time_parallel = 0.0;
    for(run = 0; run < NUM_RUNS; run++){
        //start time
        wstart=omp_get_wtime();
        #pragma simd collapse (2)
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                T[i][j] = matrix[j][i];
            }
        }
        
        //end time
        wend=omp_get_wtime();
        wtimeImp = wend - wstart;
        total_time_parallel += wtimeImp;
    }
    avg_time_imp = total_time_parallel / NUM_RUNS;
    
    
    return T;
}

#pragma GCC optimize("O2")
bool checkSymImp(double **matrix){
    int i, j;
    bool sym = true;
    
    #pragma unroll(4) collapse(2) 
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(matrix[i][j] != matrix[j][i]){
                sym = false;
            }
        }
    }
    
    return sym;
    
}

double** matTransposeOMP(double **matrix){
    int i, j;
    double **T = (double**)malloc(n*sizeof(double*));
    int num_threads;
    int run;

    for(i=0; i<n; i++){
        T[i] = (double*)malloc(n*sizeof(double)); 
    }
#ifdef _OPENMP  
    // Test performance with thread counts of 1, 2, 4, 8, 16, 32, and 64
    printf("Num_Threads | Avg_Parallel_Time | Avg_Speedup | Avg_Efficiency | Bandwidth \n"); 
    for(num_threads=2; num_threads <=OMPTHREADS; num_threads *= 2){
        
        omp_set_num_threads(num_threads);
        total_time_parallel = 0.0;
        
        for(run = 0; run < NUM_RUNS; run++){
            //start time
            wstart=omp_get_wtime();
            
            #pragma omp parallel for collapse(2) schedule (guided, 2)
            for(i=0; i<n; i++){
                for(j=0; j<n; j++){
                    T[i][j] = matrix[j][i];
                }
            }
            
            //end time
            wend=omp_get_wtime();
            wtimePar = wend - wstart;
            total_time_parallel += wtimePar;
        }
        
        // Compute average parallel time, speedup, and efficiency
        avg_time_parallel = total_time_parallel / NUM_RUNS;
        avg_speedup = avg_time_seqT / avg_time_parallel;
        avg_efficiency = avg_speedup / num_threads;
        bandwidthPar = data_transferred / (avg_time_parallel * 1e9);  // GB/s
        
        printf("%11d | %17f | %11f | %13.2f%% | %9.4f\n", 
               num_threads, avg_time_parallel, avg_speedup, avg_efficiency * 100, bandwidthPar);
        
     }
#else
    
    printf("error _OPENMP is not defined"); 
    
#endif

    return T;
}


bool checkSymOMP(double **matrix) {
    int i, j;
    bool sym = true;

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) shared(sym) schedule (guided, 2)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                sym = false;
            }   
        }
    }
    
    return sym;
#else
    printf("error _OPENMP is not defined\n");
    return false;
#endif
}
