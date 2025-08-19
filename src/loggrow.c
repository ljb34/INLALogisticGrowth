#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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

void Lmat(double growth, double carry_cap, double move_const, double step_size,
	double* linpoint, int ns, int nt, inla_cgeneric_mat_tp* CinvG, double* result) { /* result is nsnt x nsnt double matrix in inla_cgeneric_mat_tp form*/
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
		result[i * ntotal + i] = 1 / (step_size)+a_array[i];
	}

	//subdiagonal
	for (int i = 0; i < ns * (nt - 1); i++) {
		result[(ns + i) * ntotal + i] = -1 / (step_size);
	}

	//CinvG part
	for (int t = 1; t < nt; t++) {
		for (int j = 0; j < ns; j++) {
			for (int i = 0; i < ns; i++) {
				result[(ns * t + i) * ntotal + t * ns + j] += move_const * CinvG->x[i * ns + j];
			}
		}
	}
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
    double exponent = -pow(x - mean, 2.0) / (2.0 * stddev * stddev);
    return log(coefficient * exp(exponent));
}


double* inla_cgeneric_loggrow_model(inla_cgeneric_cmd_tp cmd, double* theta, inla_cgeneric_data_tp* data) {
    // this reimplement `inla.rgeneric.iid.model` using cgeneric
    double* ret = NULL, growth = (theta ? theta[0] : NAN), carry_cap = (theta ? exp(theta[1]) : NAN), move_const = (theta ? theta[2] : NAN), sigma = (theta ? exp(theta[3]) : NAN); //interpret.theta equivalent
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
    assert(!strcasecmp(data->doubles[0]->name, "step_size"));
    double step_size = data->doubles[0]->doubles[0];
    assert(step_size > 0);

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

    assert(!strcasecmp(data->mats[0]->name, "CinvG"));
    inla_cgeneric_mat_tp* CinvG = data->mats[0];
    assert(CinvG->nrow == ns);

    assert(!strcasecmp(data->mats[1]->name, "prior_variance"));
    inla_cgeneric_mat_tp* prior_variance = data->mats[1];
    assert(prior_variance->nrow == ns); 

    //different outputs for the commands supplied
    switch (cmd) { 
    case INLA_CGENERIC_VOID:
    {
        assert(!(cmd == INLA_CGENERIC_VOID));
    }
    break;
    case INLA_CGENERIC_GRAPH:
    {
        // return a vector of indices with format
        // c(N, M, ii, jj)
        // where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
        // so like the loop
        // for i=0, ...
        // for j=i, ...
        // G_ij =
        // and M is the total length while N is the dimension
        ret = Calloc(2 + N * (N + 1) / 2, double);
        assert(ret);
        ret[0] = N; /* dimension */
        ret[1] = N * (N + 1) / 2; /* number of (i <= j) */
        int idx = 2; // Start after N and M
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                ret[idx] = i; /* ii */
                ret[N*(N+1)/2 + idx] = j; /* jj */
                idx++;
            }
        }
    }
    break;
    case INLA_CGENERIC_Q:
    {
        // return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
        int M = N;
        ret = Calloc(2 + N * N, double);

        inla_cgeneric_mat_tp* L_mat = malloc(sizeof(inla_cgeneric_mat_tp));
        L_mat->x = calloc(N * N, sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        Lmat(growth, carry_cap, move_const, step_size, linpoint->doubles, ns, nt, CinvG, L_mat->x);

        //print L_mat
        /*for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < N; jj++) {
                printf("%f \t", L_mat->x[ii * N + jj]);
            }
            printf("\n");
        }*/ //checked 14/8 and it comes out correctly

        //Noise matrix
        inla_cgeneric_mat_tp* noise = malloc(sizeof(inla_cgeneric_mat_tp));
        noise->x = calloc(N * N, sizeof(double));
        noise->nrow = N;
        noise->ncol = N;
        //initial year variance
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                noise->x[i * N + j] = prior_variance->x[i * N + j];
            }
        }

        //Other years, noise on diagonal
        for (int i = ns; i < ns * nt; i++) {
            noise->x[i * N + i] = (sigma * step_size) * (sigma * step_size);
        }

        /*for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < N; jj++) {
                printf("%f \t", noise->x[ii * N + jj]);
            }
            printf("\n");
        }*/ //noise checked 15/8, comes out fine
        int* ipiv = malloc(ns * nt * sizeof(int));
        int lda = N;
        int ldb = N;
        int nrhs = N;
        int info;

        double* A = malloc(N * N * sizeof(double));
        double* B = malloc(N * N * sizeof(double));
        memcpy(A, noise->x, N * N * sizeof(double));
        // memcpy(B, L_mat->x, N * N * sizeof(double));
         // Convert L_mat->x (row-major) to B (column-major)
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                B[col * N + row] = L_mat->x[row * N + col];
            }
        }
        //Compute Noise^-1 * L
        dgesv_(&N, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        if (info != 0) {
            fprintf(stderr, "dgesv failed, info = %d\n", info);
            free(A); free(B); free(ipiv);
            free(L_mat->x); free(L_mat);
            free(noise->x); free(noise);
            return NULL;
        }
        //B is in column-major order now, print it
        /*for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < N; jj++) {
                printf("%f \t", B[jj * N + ii]);
            }
            printf("\n");
        }*/
        //Compute L^T*Noise^-1 * L
        double* out = calloc(N * N, sizeof(double));
        char transA = 'N';  // Don't Transpose A b/c in row major order
        char transB = 'N';  // No transpose B
        double one = 1, zero = 0;

        dgemm_(&transA, &transB,
            &N, &N, &N,
            &one, L_mat->x, &lda,
            B, &ldb,
            &zero, out, &N);
        free(A);
        free(B);
        free(ipiv);

        ret[0] = -1; /* REQUIRED! */
        ret[1] = N; /* number of rows */
        // Fill in ret with upper triangular part of out
        int idx = 2; // Start after -1 and M
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                ret[idx++] = out[i * N + j];
            }
        }
    }
    break;
    case INLA_CGENERIC_MU:
    {
        // return (N, mu)

        double* ret = Calloc(1 + N, double);
        assert(ret);
        ret[0] = N; /* dimension */

        inla_cgeneric_mat_tp* L_mat = malloc(sizeof(inla_cgeneric_mat_tp));
        L_mat->x = calloc(N * N, sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        Lmat(growth, carry_cap, move_const, step_size, linpoint->doubles, ns, nt, CinvG, L_mat->x);

        inla_cgeneric_vec_tp* rvector = malloc(sizeof(inla_cgeneric_vec_tp));
        rvector->doubles = calloc(N, sizeof(double));
        rvector->len = N;
        r_vector(growth, carry_cap, move_const, linpoint->doubles, mag_grad_sq->doubles, ns, nt, rvector->doubles);
        for (int i = 0; i < ns; i++) {
            rvector->doubles[i] = prior_mean->doubles[i];
        }

        //calculate L_mat^-1 * rvector
        int* ipiv = malloc(ns * nt * sizeof(int));
        int lda = N;
        int ldb = N;
        int nrhs = N;
        int info;
        double* A = malloc(N * N * sizeof(double));
        double* B = malloc(N * sizeof(double));

        //convert L_mat->x (row-major) to A (column-major)
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                A[col * N + row] = L_mat->x[row * N + col];
            }
        }

        memcpy(B, rvector->doubles, N * sizeof(double));

        dgesv_(&N, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        if (info != 0) {
            printf("dgesv failed, info = %d\n", info);
            free(A); free(B); free(ipiv);
            free(L_mat->x); free(L_mat);
            free(rvector->doubles); free(rvector);
            return NULL;
        }

        for (int i = 0; i < N; i++) {
            ret[i + 1] = B[i]; // Fill in mu
        }
        free(A);
        free(B);
        free(ipiv);
        free(L_mat->x);
        free(L_mat);
        free(rvector->doubles);
        free(rvector);
    }
    break;

    case INLA_CGENERIC_INITIAL:
    {
        // return c(P, initials)
        // where P is the number of hyperparameters
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
            ret[4] = log(5);
        }
        else {
            ret[4] = initial_sigma;
        }
    }
    break;
    case INLA_CGENERIC_LOG_PRIOR:
    {
        // return c(LOG_PRIOR)
        double* ret = Calloc(1, double);

        ret[0] = normal_pdf_log(growth, pgrowth->doubles[0], pgrowth->doubles[1]) +
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
