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
					result[(ns * t + j) * ntotal + t * ns + i] += move_const * CinvG->x[j * ns + i]; //Lmat is col-major but CinvG is row-major
                }
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

    assert(!strcasecmp(data->smats[2]->name, "C"));
    inla_cgeneric_smat_tp* C = data->smats[2];
    assert(C->nrow == ns);

    assert(!strcasecmp(data->smats[3]->name, "G"));
    inla_cgeneric_smat_tp* G = data->smats[3];
    assert(G->nrow == ns);

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
        // return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
        if (debug > 0) {
            printf("INLA_CGENERIC_Q\n");
        }
        int M = 2 * ns * ns + ns + (nt - 2) * (ns * ns + ns * (ns + 1) / 2);
        //printf("M: %d\n", M);
        ret = Calloc(2 +M, double);

        inla_cgeneric_mat_tp* L_mat = malloc(sizeof(inla_cgeneric_mat_tp));
        L_mat->x = calloc(N * N, sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        Lmat(growth, carry_cap, move_const, timestep, linpoint->doubles, ns, nt, CinvG, L_mat->x);

        int* ipiv = malloc(ns * nt * sizeof(int));
        int lda = N;
        int ldb = N;
        int nrhs = N;
        int info;

        double* B = calloc(N * N, sizeof(double));
        //Compute Noise * L
		//initial ns x ns prior_precision block - prior precision * identity so can just copy prior_precision to B
		//if sparse prior_precision
        if (prior_precision->n != ns * ns) {
            // sparse prior_precision: apply its nonzeros 
            for (int k = 0; k < prior_precision->n; k++) {
                int ii = prior_precision->i[k]; 
                int jj = prior_precision->j[k]; 
                double pv = prior_precision->x[k];
				B[jj * N + ii] = pv;
            }
        }
        else { //if dense prior precision
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
					B[j * N + i] = prior_precision->x[j * ns + i]; 
                }
            }
        }
        /*Rest of B=Q*L */
        double one = 1, zero = 0;
        //Calculate sigma**2/h * (C+gG)
        double* a_array = malloc(ns*nt * sizeof(double));
        /*a_func(growth, carry_cap, linpoint->doubles, ns, nt, a_array);
		double mean_a = 0.0;
        for (int i = 0; i < ns*nt; i++) {
            mean_a += a_array[i];
        }
        mean_a = mean_a/(ns*nt);
		free(a_array);
		double g = move_const / mean_a; */
        double g = move_const;
		double* Qblock = calloc(ns * ns, sizeof(double));
		if ((C->n == ns) & (G->n == ns)) { //if dense C and G
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    Qblock[j * ns + i] = sigma * sigma * (C->x[j * ns + i] + g * G->x[j * ns + i]) / timestep;
                }
            }
        }
        else //C and G both sparse, non zero entries don't line up
        { //calculate C+ gG
            //first copy C to Qblock
            for (int k = 0; k < C->n; k++) {
                int ii = C->i[k];
                int jj = C->j[k];
                double cv = C->x[k];
                Qblock[jj * ns + ii] = cv;
            }
            //then add gG to Qblock
            for (int k = 0; k < G->n; k++) {
                int ii = G->i[k];
                int jj = G->j[k];
                double gv = G->x[k];
                Qblock[jj * ns + ii] += g * gv;
            }
            //finally scale by sigma^2 / timestep
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    Qblock[j * ns + i] = sigma * sigma * Qblock[j * ns + i] / timestep;
                }
			}
        }
        
		//main diagonal blocks
        for (int t = 1; t < nt; t++) {

			//extract t-th block of L
			double* Lblock = malloc(ns * ns * sizeof(double));
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    Lblock[j * ns + i] = L_mat->x[(t * ns + j) * N + t * ns + i];
                }
            }
			//multiply Qblock * Lblock and store in B
			double* temp = malloc(ns * ns * sizeof(double));
			char transA = 'N';
			char transB = 'N';
			dgemm_(&transA, &transB, &ns, &ns, &ns, &one, Qblock, &ns, Lblock, &ns, &zero, temp, &ns);
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    B[(t * ns + j) * N + t * ns + i] = temp[j * ns + i];
                }
            }

			//clear memory
			free(Lblock);
			free(temp);
        }

		//sub diagonal blocks -1/timestep*Qblock
        for (int t = 1; t < nt - 1; t++) {
            for(int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    B[(t * ns + j) * N + (t + 1) * ns + i] = -1.0 / timestep * Qblock[j * ns + i];
                }
			}
        }

		free(Qblock);

        /* Now compute out = L_mat^T * B */
        double* out = calloc(N * N, sizeof(double));
        char transL = 'T';
		char transB = 'N';
        dgemm_(&transL, &transB, &N, &N, &N, &one, L_mat->x, &lda, B, &ldb, &zero, out, &N);

        free(L_mat->x);
        free(L_mat);
        free(B);
        free(ipiv);

        ret[0] = -1; /* REQUIRED! */
        ret[1] = M; 

		//fill in ret with non zero parts of out 


		int idx = 2; // Start after -1 and M
		//first year only has two blocks
        for (int i = 0; i < ns; i++) {
            for (int j = i; j < 2*ns; j++) {
                ret[idx++] = out[j * N + i];
            }
        }
		//middle years have three blocks
        for (int k = 1; k < nt - 1; k++) {
            for (int i = k * ns; i < (k + 1) * ns; i++) {
                for (int j = i; j < (k + 2) * ns; j++) {
                    ret[idx++] = out[j * N + i];
                }
            }
        }

        //final year has two blocks
        
        for (int i = (nt - 1) * ns; i < nt * ns; i++) {
            for (int j = i; j < nt * ns; j++) {
                ret[idx++] = out[j * N + i];
            }
        }
		free(out);

        if (idx - 2 != M) {
            fprintf(stderr, "Q filled %d values, expected %d\n", idx - 2, M);
            abort();
        }
    }
    break;
    case INLA_CGENERIC_MU:
    {
        // return (N, mu)
        if (debug > 0) {
            printf("INLA_CGENERIC_MU\n");
        }
        ret = Calloc(1 + N, double);
        assert(ret);
        ret[0] = N; /* dimension */

        inla_cgeneric_mat_tp* L_mat = malloc(sizeof(inla_cgeneric_mat_tp));
        L_mat->x = calloc(N * N, sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        Lmat(growth, carry_cap, move_const, timestep, linpoint->doubles, ns, nt, CinvG, L_mat->x);

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
        int nrhs = 1;
        int info;
        double* A = malloc(N * N * sizeof(double));
        double* B = malloc(N * sizeof(double));
		memcpy(A, L_mat->x, N* N * sizeof(double));
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
