#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "cgeneric.h"
#define Calloc(n_, type_) (type_ *)calloc((n_), sizeof(type_))
#if !defined(iszero)
#ifdef __SUPPORT_SNAN__
#define iszero(x) (fpclassify(x) == FP_ZERO)
#else
#define iszero(x) (((__typeof(x))(x)) == 0)
#endif
#endif

extern void dgesv_(int* n, int* nrhs, double* A, int* lda,
	int* ipiv, double* B, int* ldb, int* info);

extern void dgemm_(const char* TRANSA, const char* TRANSB,
	const int* M, const int* N, const int* K,
	const double* alpha,
	const double* A, const int* lda,
	const double* B, const int* ldb,
	const double* beta,
	double* C, const int* ldc);

void a_func(double growth, double carry_cap,
	double* linpoint, int ns, int nt, double* result) { /*important! must give pointer for where to export result*/
	for (int i = 0; i < ns * nt; i++) {
		result[i] = growth * exp(linpoint[i]) / carry_cap;
	}
}
//Full dense version of Lmat, output is in COLUMN MAJOR
void Lmat(double growth, double carry_cap, double move_const, double timestep,
	double* linpoint, int ns, int nt, inla_cgeneric_smat_tp* CinvG, double* result) { /* result is nsnt x nsnt double matrix in inla_cgeneric_mat_tp form*/
	//identity sub matrix in first block
	int ntotal = ns * nt;
	for (int i = 0; i < ns; i++) {
		result[i * ntotal + i] = 1;
	}

	//main diagonal
	double* a_array = malloc(ntotal * sizeof(double));
	a_func(growth, carry_cap,
		linpoint, ns, nt, a_array);
	for (int i = ns; i < ntotal; i++) {
		result[i * ntotal + i] = a_array[i] + (1 /timestep);
	}

	//subdiagonal
	for (int i = 0; i < ns * (nt - 1); i++) {
		result[i * ntotal + i + ns] = -1 / timestep;
	}

	//CinvG part
    if(CinvG->n != ns * ns) {
        for(int idx = 0; idx < CinvG->n; idx++) {
            int i = CinvG->i[idx];
            int j = CinvG->j[idx];
            for (int t = 1; t < nt; t++) {
                result[(ns * t + j) * ntotal + t * ns + i] += move_const * CinvG->x[idx];
            }
		}
    }
    else {
        for (int t = 1; t < nt; t++) {
            for (int j = 0; j < ns; j++) {
                for (int i = 0; i < ns; i++) {
					result[(ns * t + j) * ntotal + t * ns + i] += move_const * CinvG->x[i * ns + j]; //Lmat is col-major but CinvG is row-major
                }
            }
        }
    }
	free(a_array);

}

// Sparse version of Lmat, output is in COLUMN MAJOR
void Lmat_sparse(double growth, double carry_cap, double move_const, double timestep,
    double* linpoint, int ns, int nt, inla_cgeneric_smat_tp* CinvG, inla_cgeneric_smat_tp* result) {
	//identity sub matrix in first block
    for (int i = 0; i < ns; i++) {
		result->i[i] = i;
		result->j[i] = i;
		result->x[i] = 1;
    }
    
	//subdiagonal
    for(int i = 0; i < ns * (nt - 1); i++) {
        result->i[ns * nt + i] = i;
        result->j[ns * nt + i] = i + ns;
        result->x[ns * nt + i] = -1 / timestep;
	}
	//Main diagonal block - CinvG + diag(a_array + 1/timestep)
	double* a_array = malloc(ns * nt * sizeof(double));
    a_func(growth, carry_cap,
        linpoint, ns, nt, a_array);
           
	int offset = ns * nt + ns * (nt - 1);
    if (CinvG->n != ns * ns) {
		printf("Sparse CinvG not supported in Lmat_sparse yet\n");
    }
    for (int t = 1; t < nt; t++) {
        for (int i = t * ns; i < (t + 1) * ns; i++) {
            for (int j = t * ns; j < (t + 1) * ns; j++) {
                if (i == j) { //on diagonal, include a_array + 1/timestep
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = move_const * CinvG->x[(i - t * ns) * ns + (j - t * ns)] + a_array[i] + (1 / timestep);
                }
                else {
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = move_const * CinvG->x[(i - t * ns) * ns + (j - t * ns)];
				}
				offset++;
            }
        }
    }
	free(a_array);

}

// Block version of Lmat, output is in COLUMN MAJOR
void Lmat_block(double growth, double carry_cap, double move_const, double timestep,
    double* linpoint, int ns, int nt, inla_cgeneric_smat_tp* CinvG, inla_cgeneric_smat_tp* result) {
    //identity sub matrix in first block
	int offset = 0;
    for (int i = 0; i < ns; i++) { //identity sub matrix in first block
        for (int j = 0; j < ns; j++) {
            if (i == j) {
                result->i[offset] = i;
                result->j[offset] = j;
                result->x[offset] = 1;
            }
            else {
                result->i[offset] = i;
                result->j[offset] = i;
                result->x[offset] = 0;
            }
            offset++;
        }
    }

    //subdiagonal
    for (int t = 1; t < nt - 1; t++) {
        for (int i = (t - 1) * ns; i < t * ns; i++) {
            for (int j = t * ns; j < (t + 1) * ns; j++) {
                if (i + ns == j) {
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = -1 / timestep;
                }
                else {
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = 0;
                }
                offset++;
            }
        }
    }
    //Main diagonal block - CinvG + diag(a_array + 1/timestep)
    double* a_array = malloc(ns * nt * sizeof(double));
    a_func(growth, carry_cap,
        linpoint, ns, nt, a_array);

    if (CinvG->n != ns * ns) {
        printf("Sparse CinvG not supported in Lmat_sparse yet\n");
    }
    for (int t = 1; t < nt; t++) {
        for (int i = t * ns; i < (t + 1) * ns; i++) {
            for (int j = t * ns; j < (t + 1) * ns; j++) {
                if (i == j) { //on diagonal, include a_array + 1/timestep
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = move_const * CinvG->x[(i - t * ns) * ns + (j - t * ns)] + a_array[i] + (1 / timestep);
                }
                else {
                    result->i[offset] = i;
                    result->j[offset] = j;
                    result->x[offset] = move_const * CinvG->x[(i - t * ns) * ns + (j - t * ns)];
                }
                offset++;
            }
        }
    }
    free(a_array);
   
}
void r_vector(double growth, double carry_cap, double move_const,
    double* linpoint, double* mag_grad_sq, int ns, int nt, double* result) {
    for (int i = 0; i < ns; i++) {
        for (int t = 0; t < nt; t++) {
            int idx = t * ns + i;
            //mag_grad_sq = grad[2 * idx] * grad[2 * idx] + grad[2 * idx + 1] * grad[2 * idx + 1];
            double lp = linpoint[idx];
            result[idx] = growth * exp(lp) * (lp - 1) / carry_cap + growth - move_const * mag_grad_sq[idx];
        }
    }
}

double normal_pdf_log(double x, double mean, double stddev) {
    double coefficient = 1.0 / (stddev * sqrt(2.0 * M_PI));
    double exponent = -(x - mean)*(x - mean)/ (2 * stddev * stddev);
    return log(coefficient * exp(exponent));
}


double* inla_cgeneric_loggrow_model(inla_cgeneric_cmd_tp cmd, double* theta, inla_cgeneric_data_tp* data) {
    // this reimplement `inla.rgeneric.iid.model` using cgeneric
    double* ret = NULL, growth = (theta ? exp(theta[0]) : NAN), carry_cap = (theta ? exp(theta[1]) : NAN), move_const = (theta ? theta[2] : NAN), sigma = (theta ? exp(theta[3]) : NAN); //interpret.theta equivalent
    assert(!strcasecmp(data->ints[0]->name, "n")); // this will always be the case
    int N = data->ints[0]->ints[0]; // this will always be the case
    assert(N > 0);

    // Read parameters from data
    assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
    int debug = data->ints[1]->ints[0];
    if (debug > 0) debug = 1;

    assert(!strcasecmp(data->ints[2]->name, "ns"));
    int ns = data->ints[2]->ints[0];
    assert(ns > 0);

    assert(!strcasecmp(data->ints[3]->name, "nt"));
    int nt = data->ints[3]->ints[0];
    assert(nt > 0);

    //Pre calculated information
    assert(!strcasecmp(data->doubles[0]->name, "timestep"));
    double timestep = data->doubles[0]->doubles[0];
    assert(timestep > 0);

    assert(!strcasecmp(data->doubles[1]->name, "linpoint"));
    inla_cgeneric_vec_tp* linpoint = data->doubles[1];
    assert(linpoint->len == ns * nt);

    /*assert(!strcasecmp(data->doubles[2]->name, "CinvG"));
    inla_cgeneric_mat_tp* CinvG = data->doubles[2];
    assert(CinvG->nrow == ns);

    assert(!strcasecmp(data->doubles[3]->name, "prior_variance"));
    inla_cgeneric_mat_tp* prior_variance = data->doubles[3];
    assert(prior_variance->nrow == ns);*/ //matrices not doubles

    assert(!strcasecmp(data->doubles[2]->name, "mag_grad_sq"));
    inla_cgeneric_vec_tp* mag_grad_sq = data->doubles[2];
    assert(mag_grad_sq->len == ns * nt);

    assert(!strcasecmp(data->doubles[3]->name, "prior_mean"));
    inla_cgeneric_vec_tp* prior_mean = data->doubles[3];
    assert(prior_mean->len == ns);

    //initial values
    assert(!strcasecmp(data->doubles[4]->name, "initial_growth"));
    double initial_growth = data->doubles[4]->doubles[0];

    assert(!strcasecmp(data->doubles[5]->name, "initial_carry_cap"));
    double initial_carry_cap = data->doubles[5]->doubles[0];

    assert(!strcasecmp(data->doubles[6]->name, "initial_move_const"));
    double initial_move_const = data->doubles[6]->doubles[0];

    assert(!strcasecmp(data->doubles[7]->name, "initial_sigma"));
    double initial_sigma = data->doubles[7]->doubles[0];



    //prior paramters
    assert(!strcasecmp(data->doubles[8]->name, "pgrowth"));
    inla_cgeneric_vec_tp* pgrowth = data->doubles[8];
    assert(pgrowth->len == 2);

    assert(!strcasecmp(data->doubles[9]->name, "pcc"));
    inla_cgeneric_vec_tp* pcc = data->doubles[9];
    assert(pcc->len == 2);

    assert(!strcasecmp(data->doubles[10]->name, "pmove"));
    inla_cgeneric_vec_tp* pmove = data->doubles[10];
    assert(pmove->len == 2);

    assert(!strcasecmp(data->doubles[11]->name, "psigma"));
    inla_cgeneric_vec_tp* psigma = data->doubles[11];
    assert(psigma->len == 2);

	//matrices
    assert(!strcasecmp(data->smats[0]->name, "CinvG"));
    inla_cgeneric_smat_tp* CinvG = data->smats[0];
    assert(CinvG->nrow == ns);

    assert(!strcasecmp(data->smats[1]->name, "prior_precision"));
    inla_cgeneric_smat_tp* prior_precision = data->smats[1];
    assert(prior_precision->nrow == ns); 

    //different outputs for the commands supplied
    switch (cmd) { 
    case INLA_CGENERIC_VOID:
    {
        assert(!(cmd == INLA_CGENERIC_VOID));
    }
    break;
    case INLA_CGENERIC_GRAPH:
    {   
		printf("Calculating GRAPH \n");
        // return a vector of indices with format
        // c(N, M, ii, jj)
        // where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
        // so like the loop
        // for i=0, ...
        // for j=i, ...
        // G_ij =
        // and M is the total length while N is the dimension
        int M =  2* ns * ns + ns  +(nt-2)*(ns * ns + ns * (ns+1)/2);
		//printf("M: %d\n", M);
        ret = calloc(2 + 2*M, sizeof(double));
        assert(ret);
        ret[0] = N; /* dimension */
        ret[1] = M; /* number of (i <= j) */
        int idx = 2; // Start after N and M
        //first year only has two blocks
        for (int i = 0; i < ns;i++) {
            for (int j = i; j < 2*ns; j++) {
				ret[idx] = i; /* ii */
				ret[M + idx] = j; /* jj */
                idx++;
            }
        }
		//middle years have three blocks
        for (int k = 1; k < nt-1; k++) {
            for(int i = k*ns; i < (k+1)*ns; i++) {
                for(int j = i; j < (k+2)*ns; j++) {
                    ret[idx] = i; /* ii */
                    ret[M + idx] = j; /* jj */
					idx++;
                }
			}
        }
		//final year only has two blocks
        for (int i = (nt - 1) * ns; i < nt * ns; i++) {
            for (int j = i; j < nt * ns; j++) {
                ret[idx] = i; /* ii */
                ret[M + idx] = j; /* jj */
                idx++;
            }
        }
        //for (int i = 0; i < N; i++) {
        //    for (int j = i; j < N; j++) {
        //        ret[idx] = i; /* ii */
        //        ret[N*(N+1)/2 + idx] = j; /* jj */
        //        idx++;
        //    }
        //}
        if (idx - 2 != M) {
            fprintf(stderr, "GRAPH produced %d pairs, expected %d\n", idx - 2, M);
            abort();
        }

    }
    break;
    case INLA_CGENERIC_Q:
    {
        printf("Calculating Q \n");
        // return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
        if (debug > 0) {
            printf("INLA_CGENERIC_Q\n");
        }
        int M = 2 * ns * ns + ns + (nt - 2) * (ns * ns + ns * (ns + 1) / 2);
        //printf("M: %d\n", M);
        ret = Calloc(2 +M, double);

        inla_cgeneric_smat_tp* L_mat = malloc(sizeof(inla_cgeneric_smat_tp));
        L_mat->x = malloc(ns*ns + 2*ns*ns*(nt-1)*sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        L_mat->i = malloc(L_mat->n * sizeof(int));
        L_mat->j = malloc(L_mat->n * sizeof(int));
		L_mat->n = ns * ns + 2 * ns * ns * (nt - 1); //number of nonzeros in L
        Lmat_block(growth, carry_cap, move_const, timestep, linpoint->doubles, ns, nt, CinvG, L_mat);

		//Make copy of L_mat
		
        inla_cgeneric_smat_tp* B = malloc(sizeof(inla_cgeneric_smat_tp));
        B->nrow = L_mat->nrow;
        B->n = L_mat->n;
        B->i = malloc(L_mat->n * sizeof(int));
        B->j = malloc(L_mat->n * sizeof(int));
        B->x = malloc(L_mat->n * sizeof(double));
        memcpy(B->i, L_mat->i, L_mat->n * sizeof(int));
        memcpy(B->j, L_mat->j, L_mat->n * sizeof(int));
        memcpy(B->x, L_mat->x, L_mat->n * sizeof(double));
  
        //Compute Noise * L
		//scale everything but first ns rows by timestep / (sigma^2)
		double scale = timestep / (sigma * sigma);
            for (int i = ns*ns; i < B->n; i++) {
				B->x[i] *= scale;
            }

		//initial ns x ns prior_precision block - prior precision * identity so can just copy prior_precision to B
		//if sparse prior_precision
            if (prior_precision->n != ns * ns) {
                // sparse prior_precision: apply its nonzeros 
                for (int k = 0; k < prior_precision->n; k++) {
                    int ii = prior_precision->i[k];
                    int jj = prior_precision->j[k];
                    //find corresponding entry in B
                    int found = 0;
                    for (int idx = 0; idx < ns * ns; idx++) {
                        if (B->i[idx] == jj && B->j[idx] == ii) { //B is col-major but prior_precision is row-major
                            B->x[idx] += prior_precision->x[k];
                            found = 1;
                            break;
                        }
                    }
                    if (found != 1) {
                        fprintf(stderr, "Could not find matching entry in B for prior precision at (%d, %d)\n", ii, jj);
                        abort();
                    }
                }
            }
        else {
            // dense prior_precision
            for (int i = 0; i < ns; i++) {
				for (int j = 0; j < ns; j++) {
					B->x[i * ns + j] += prior_precision->x[j * ns + i]; //B is col-major but prior_precision is row-major
                    }
            }
		}

        
        ret[0] = -1; /* REQUIRED! */
        ret[1] = M;
        int idx = 2; // Start after -1 and M
		// Compute out = L_mat^T * B block by block

        //first ns rows
        for (int i = 0; i < ns; i++) {
            //first column block
            for (int j = i; j < ns; j++) { //ret needs to be upper triangular
                if (i == j) {
                    ret[idx] = B->x[j * ns + i] + 1 / timestep * sigma * sigma;

                }
                else {
                    ret[idx] = B->x[j * ns + i];
                }
                idx++;
            }
            //second column block
            for (int j = i; j < 2 * ns; j++) {
                ret[idx] = -1 / timestep * B->x[j * ns + i];
                idx++;
            }
        }

		//middle blocks
		//set up for dgemm
        const char transL = 'T';   // L^T
        const char transB = 'N';   // B
		const int nrowsC = ns;     //size of output block 
        const int lda = ns;         // leading dimension of L
        const int ldb = ns;         // leading dimension of B
        const double alpha = 1.0;
        const double beta = 0.0;

		int diag_start = ns * ns + (nt - 1) * ns * ns; //start of diagonal blocks in B
        for (int t = 1; t < nt - 2; t++) {
            //diagonal
            // extract L  & B blocks
            double* L_block = malloc(ns * ns * sizeof(double));
            double* B_block = malloc(ns * ns * sizeof(double));
            for (int i = t * ns; i < (t + 1) * ns; i++) {
                for (int j = t * ns; j < (t + 1) * ns; j++) {
                    L_block[(j - t * ns) * ns + (i - t * ns)] = L_mat->x[diag_start + (t - 1) * ns * ns + (j - t * ns) * ns + (i - t * ns)]; //t-th block row, t-th block column
                    B_block[(j - t * ns) * ns + (i - t * ns)] = B->x[diag_start + (t - 1) * ns * ns + (j - t * ns)*ns + (i - t * ns)];
                }
            }
            double* C_block = malloc(ns * ns * sizeof(double));
            dgemm_(&transL, &transB, &nrowsC, &nrowsC, &nrowsC,
                &alpha, L_block, &lda,
                B_block, &ldb,
				&beta, C_block, &nrowsC);

            //fill in the t-th row of ret
            for(int i = t * ns; i < (t + 1) * ns; i++){
                for(int j = i; j < (t + 1) * ns; j++){
					ret[idx] = C_block[(j - t * ns) * ns + (i - t * ns)] + 1/timestep*sigma*sigma;
                    idx++;
				}
				//off diagonal block
                for(int j = (t + 1) * ns; j < (t + 2) * ns; j++){
					ret[idx] = ( - 1 / timestep) * B->x[diag_start + t * ns * ns + (j - (t + 1) * ns) * ns + (i - t * ns)]; //get the t+1,t+1 block from B
					idx++;
				}
             }
            free(L_block);
            free(B_block);
			free(C_block);
        }

		//final block - diagonal only
        double* L_block = malloc(ns * ns * sizeof(double));
        double* B_block = malloc(ns * ns * sizeof(double));
        for (int i = (nt - 1) * ns; i < nt * ns; i++) {
            for (int j = (nt - 1) * ns; j < nt * ns; j++) {
                L_block[(j - (nt - 1) * ns) * ns + (i - (nt - 1) * ns)] = L_mat->x[diag_start + (nt - 2) * ns * ns + (j - (nt-1) * ns) * ns + (i - (nt-1) * ns)];
                B_block[(j - (nt - 1) * ns) * ns + (i - (nt - 1) * ns)] = B->x[diag_start + (nt - 2) * ns * ns + (j - (nt - 1) * ns) * ns + (i - (nt - 1) * ns)];
            }
        }
            double* C_block = malloc(ns * ns * sizeof(double));
            dgemm_(&transL, &transB, &nrowsC, &nrowsC, &nrowsC,
                &alpha, L_block, &lda,
                B_block, &ldb,
                &beta, C_block, &nrowsC);
			//fill in the t-th row of ret
            for(int i = (nt - 1) * ns; i < nt * ns; i++){
                for(int j = i; j < nt * ns; j++){
                    ret[idx] = C_block[(j - (nt - 1) * ns) * ns + (i - (nt - 1) * ns)] + 1/timestep*sigma*sigma;
                    idx++;
                }
             }
        free(L_block);
		free(B_block);
		free(C_block);
        free(L_mat->x);
        free(L_mat->i);
        free(L_mat->j);
        free(L_mat);

        free(B->x);
		free(B->i);
		free(B->j);
        free(B);
		
        if (idx - 2 != M) {
            fprintf(stderr, "Q filled %d values, expected %d\n", idx - 2, M);
            abort();
        }
    }
    break;
    case INLA_CGENERIC_MU:
    {   
		printf("Calculating MU \n");
        // return (N, mu)
        if (debug > 0) {
            printf("INLA_CGENERIC_MU\n");
        }
        ret = Calloc(1 + N, double);
        assert(ret);
        ret[0] = N; /* dimension */
		int idx = 1;
        inla_cgeneric_smat_tp* L_mat = malloc(sizeof(inla_cgeneric_smat_tp));
        L_mat->x = malloc(ns * ns + 2 * ns * ns * (nt - 1)*sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        L_mat->i = malloc(L_mat->n * sizeof(int));
        L_mat->j = malloc(L_mat->n * sizeof(int));
        L_mat->n = ns * ns + 2 * ns * ns * (nt - 1); //number of nonzeros in L
        Lmat_block(growth, carry_cap, move_const, timestep, linpoint->doubles, ns, nt, CinvG, L_mat);

        inla_cgeneric_vec_tp* rvector = malloc(sizeof(inla_cgeneric_vec_tp));
        rvector->doubles = calloc(N, sizeof(double));
        rvector->len = N;
        r_vector(growth, carry_cap, move_const, linpoint->doubles, mag_grad_sq->doubles, ns, nt, rvector->doubles);
        for (int i = 0; i < ns; i++) {
			ret[idx] = prior_mean->doubles[i];
            idx++;
        }

        //calculate L_mat^-1 * rvector block by block
        int* ipiv = malloc(ns * nt * sizeof(int));
        int lda = ns;
        int ldb = ns;
        int nrhs = 1;
        int info;
        int diag_start = ns * ns + (nt - 1) * ns * ns; //start of diagonal blocks in B

        for (int t = 1; t < nt; t++) {
			//extract L block
			double* L_block = malloc(ns * ns * sizeof(double));
			double* rvector_block = malloc(ns * sizeof(double));
            for (int i = t * ns; i < (t + 1) * ns; i++) {
                for (int j = t * ns; j < (t + 1) * ns; j++) {
                    L_block[(j - t * ns) * ns + (i - t * ns)] = L_mat->x[diag_start + (t - 1) * ns * ns + (j - t * ns) * ns + (i - t * ns)]; //t-th block row, t-th block column
                }
				rvector_block[i - t * ns] = rvector->doubles[i];
            }
			//solve L_block * x = rvector_block
            dgesv_(&ns, &nrhs, L_block, &lda,
				ipiv, rvector_block, &ldb, &info);
            if (info != 0) {
                fprintf(stderr, "Error in dgesv_: info = %d\n", info);
                abort();
			}
            //put into ret
            for (int i = 0; i < ns; i++) {
                ret[idx] = rvector_block[i];
                idx++;
            }
			free(L_block);
			free(rvector_block);
        }
        free(ipiv);
        free(L_mat->x);
        free(L_mat->i);
        free(L_mat->j);
        free(L_mat);
        free(rvector->doubles);
        free(rvector);
    }
    break;

    case INLA_CGENERIC_INITIAL:
    {
        // return c(P, initials)
        // where P is the number of hyperparameters
        if(debug>0) {
            printf("INLA_CGENERIC_INITIAL\n");
		}
        ret = Calloc(5, double);
        ret[0] = 4;

        if (iszero(initial_growth)) {
            ret[1] = 1;
        }
        else {
            ret[1] = initial_growth;
        }

        if (iszero(initial_carry_cap)) {
            ret[2] = log(100);
        }
        else {
            ret[2] = initial_carry_cap;
        }

        if (iszero(initial_move_const)) {
            ret[3] = 1;
        }
        else {
            ret[3] = initial_move_const;
        }

        if (iszero(initial_sigma)) {
            ret[4] = log(1);
        }
        else {
            ret[4] = initial_sigma;
        }
    }
    break;
    case INLA_CGENERIC_LOG_PRIOR:
    {
        // return c(LOG_PRIOR)
        if (debug > 0) {
            printf("INLA_CGENERIC_LOG_PRIOR\n");
        }
        ret = Calloc(1, double);

        ret[0] = normal_pdf_log(theta[0], pgrowth->doubles[0], pgrowth->doubles[1]) +
            normal_pdf_log(theta[1], pcc->doubles[0], pcc->doubles[1]) + 
            normal_pdf_log(move_const, pmove->doubles[0], pmove->doubles[1]) + 
            normal_pdf_log(theta[3], psigma->doubles[0],psigma->doubles[1]);
    }
    break;

    case INLA_CGENERIC_LOG_NORM_CONST:
    case INLA_CGENERIC_QUIT:
    default:
        break;
    }
    return(ret);
}
