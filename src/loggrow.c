#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "cgeneric.h"
#include "utils.h"
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
					result[(ns * t + j) * ntotal + t * ns + i] += move_const * CinvG->x[j * ns + i]; 
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

double sparse_get(
    int row,
    int col,
    int n_entries,
    int* I,
    int* J,
    double* X)
{
    for (int k = 0; k < n_entries; k++) {
        if (I[k] == row && J[k] == col)
            return X[k];

        // if symmetric storage:
        if (I[k] == col && J[k] == row)
            return X[k];
    }
    return 0.0;
}


double* inla_cgeneric_loggrow_model(inla_cgeneric_cmd_tp cmd, double* theta, inla_cgeneric_data_tp* data) {
    // this reimplement `inla.rgeneric.iid.model` using cgeneric
    double* ret = NULL, growth = (theta ? exp(theta[0]) : NAN), carry_cap = (theta ? exp(theta[1]) : NAN), move_const = (theta ? exp(theta[2]) : NAN), sigma = (theta ? exp(theta[3]) : NAN); //interpret.theta equivalent
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

    assert(!strcasecmp(data->ints[4]->name, "Pn"));
    int Pn = data->ints[4]->ints[0];
    assert(Pn > 0);

    assert(!strcasecmp(data->ints[5]->name, "offdn"));
	int offdn = data->ints[5]->ints[0];
	assert(offdn > 0);

	assert(!strcasecmp(data->ints[6]->name, "diagn"));
	int diagn = data->ints[6]->ints[0];
	assert(diagn > 0);

    //Non-zero locations
    assert(!strcasecmp(data->ints[7]->name, "Pi"));
    inla_cgeneric_vec_tp* Pi = data->ints[7];
    assert(Pi->len == Pn);

    assert(!strcasecmp(data->ints[8]->name, "Pj"));
    inla_cgeneric_vec_tp* Pj = data->ints[8];
    assert(Pj->len == Pn);

    assert(!strcasecmp(data->ints[9]->name, "offdi"));
    inla_cgeneric_vec_tp* offdi = data->ints[9];
    assert(offdi->len == offdn);

    assert(!strcasecmp(data->ints[10]->name, "offdj"));
    inla_cgeneric_vec_tp* offdj = data->ints[10];
    assert(offdj->len == offdn);

    assert(!strcasecmp(data->ints[11]->name, "diagi"));
    inla_cgeneric_vec_tp* diagi = data->ints[11];
    assert(diagi->len == diagn);

    assert(!strcasecmp(data->ints[12]->name, "diagj"));
    inla_cgeneric_vec_tp* diagj = data->ints[12];
    assert(diagj->len == diagn);

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
        int prev_i = -1;
        int prev_j = -1;
        // return a vector of indices with format
        // c(N, M, ii, jj)
        // where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
        // so like the loop
        // for i=0, ...
        // for j=i, ...
        // G_ij =
        // and M is the total length while N is the dimension
		int M = Pn + (nt - 1) * offdn + (nt - 1) * diagn; //total number of non zero entries in the precision matrix
        ret = calloc(2 + 2*M, sizeof(double));
        assert(ret);
        ret[0] = N; /* dimension */
        ret[1] = M; /* number of (i <= j) */
        int idx = 2; // Start after N and M
        //first year only has two blocks
        for (int i = 0; i < ns;i++) {
			//for j in Pj[Pi == i]             
            for (int k = 0; k < Pn; k++) {
                if (Pi->ints[k] == i) {
                    int j = Pj->ints[k];
                    if (j >= i) { // only include upper triangle
                        ret[idx] = i; 
                        ret[M + idx] = j; 
                        idx++;
                        if (i < prev_i || (i == prev_i && j < prev_j)) {
                            printf("GRAPH ORDER VIOLATION: (%d,%d) after (%d,%d)\n",
                                i, j, prev_i, prev_j);
                        }

                        prev_i = i;
                        prev_j = j;
                    }
                }
            }
				//for j in offdj[offdi == i], need full matrix not just upper triangle
                for (int k = 0; k < offdn; k++) {
                    if (offdi->ints[k] == i) {
                        int j = offdj->ints[k];
                        ret[idx] = i; 
                        ret[M + idx] = j + ns; 
                        idx++;
                        if (i < prev_i || (i == prev_i && j + ns < prev_j)) {
                            printf("GRAPH ORDER VIOLATION: (%d,%d) after (%d,%d)\n",
                                i, j, prev_i, prev_j);
                        }

                        prev_i = i;
                        prev_j = j + ns;
                    }
                }
        }
		//middle years have three blocks
                for (int t = 1; t < nt - 1; t++) {
                    for (int i = t * ns; i < (t + 1) * ns; i++) {
                        for(int k = 0; k < diagn; k++){
                            if (diagi->ints[k] == i - t * ns) {
                                int j = diagj->ints[k];
                                if (j >= i - t * ns) { // only include upper triangle
                                    ret[idx] = i; /* ii */
                                    ret[M + idx] = j + t*ns; /* jj */
                                    idx++;
                                    if (i < prev_i || (i == prev_i && j + t * ns < prev_j)) {
                                        printf("GRAPH ORDER VIOLATION: (%d,%d) after (%d,%d)\n",
                                            i, j, prev_i, prev_j);
                                    }

                                    prev_i = i;
                                    prev_j = j + t * ns;
                                }
                            }
						}
                        for (int k = 0; k < offdn; k++) {
                            if (offdi->ints[k] == i - t * ns) {
                                int j = offdj->ints[k];
                                ret[idx] = i; /* ii */
                                ret[M + idx] = j + (t+1)*ns; /* jj */
                                idx++;
                                if (i < prev_i || (i == prev_i && j + (t+1) * ns < prev_j)) {
                                    printf("GRAPH ORDER VIOLATION: (%d,%d) after (%d,%d)\n",
                                        i, j, prev_i, prev_j);
                                }

                                prev_i = i;
                                prev_j = j + (t+1) * ns;
                            }
                        }
                    }
                }

		//final year only has two blocks
        for (int i = (nt - 1) * ns; i < nt * ns; i++) {
            for(int k = 0; k < diagn; k++){
                if (diagi->ints[k] == i - (nt - 1) * ns) {
                    int j = diagj->ints[k];
                    if (j >= i - (nt - 1) * ns) { // only include upper triangle
                        ret[idx] = i; /* ii */
                        ret[M + idx] = j + (nt-1)*ns; /* jj */
                        idx++;
                        if (i < prev_i || (i == prev_i && j + (nt - 1) * ns < prev_j)) {
                            printf("GRAPH ORDER VIOLATION: (%d,%d) after (%d,%d)\n",
                                i, j, prev_i, prev_j);
                        }

                        prev_i = i;
                        prev_j = j + (nt - 1) * ns;
                    }
                }
			}
        }
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
            printf("nrow = %f\n", N);
        }
        int M = Pn + (nt - 1) * offdn + (nt - 1) * diagn;
        if (debug > 0) printf("M: %d\n", M);
        ret = Calloc(2 +M, double);
        ret[0] = -1; /* REQUIRED! */
        ret[1] = M;

        double g = move_const;
        //printf("move_const: %f\n", move_const);
        double* Qblock = calloc(ns * ns, sizeof(double));
        if ((C->n == ns*ns) && (G->n == ns*ns)) { //if dense C and G
            for (int i = 0; i < ns; i++) {
                for (int j = 0; j < ns; j++) {
                    Qblock[j * ns + i] = C->x[j * ns + i] + g * G->x[j * ns + i];
                    if(isnan(Qblock[j * ns + i]) & debug > 0) {
                        printf("Warning: NaN in Qblock at (%d, %d) from C and G values %f and %f\n", i, j, C->x[j * ns + i], G->x[j * ns + i]);
					}
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
                if (isnan(cv) & debug > 0) {
                    printf("Warning: NaN in C at (%d, %d)\n", ii, jj);
                }
                Qblock[jj * ns + ii] = cv;
                
            }
            //then add gG to Qblock
            for (int k = 0; k < G->n; k++) {
                int ii = G->i[k];
                int jj = G->j[k];
                double gv = G->x[k];
                if (isnan(gv) & debug > 0) {
                    printf("Warning: NaN in G at (%d, %d)\n", ii, jj);
				}
                Qblock[jj * ns + ii] += g * gv;
                
            } 
        }
        
        double one = 1, zero = 0;
        double* a_array = malloc(ns*nt * sizeof(double));
        a_func(growth, carry_cap,
            linpoint->doubles, ns, nt, a_array);
      
        //copy CinvG to new matrix and add -1/timestep to diagonal
		double* fT = calloc(ns * ns, sizeof(double));
        
        for (int k = 0; k < CinvG->n; k++) {
            int i = CinvG->i[k];
            int j = CinvG->j[k];
            double v = CinvG->x[k];

            fT[j * ns + i] = move_const * v;
            
        }
        for (int i = 0; i < ns; i++) {
            fT[i * ns + i] += a_array[ns + i] + 1.0 / timestep;
        }
        
        //calculate Q*fT and store in QfT
        double* QfT = calloc(ns * ns, sizeof(double));
        char transA = 'N';
		char transB = 'N';
		if (debug > 0) printf("dgemm step");
        dgemm_(&transA, &transB,
            &ns, &ns, &ns,
            &one,
            Qblock, &ns,
            fT, &ns,
            &zero,
            QfT, &ns);
		//start filling in ret in order of GRAPH
        int idx = 2;
        double val = 0.0;
        //first year only has two blocks
        for(int i = 0; i < ns; i++) {
            //need to extract prior_precision->x[prior_precision->i == i] for all j, then QfT[j*ns + i] for all j, then i++
            for (int k = 0; k < Pn; k++) {
                if(Pi->ints[k] == i) {
                    int j = Pj->ints[k];
                    if (j >= i) { // only include upper triangle
                        val = sparse_get(i, j,
                            prior_precision->n,
                            prior_precision->i,
                            prior_precision->j,
                            prior_precision->x) + (1 / (sigma * sigma * timestep * timestep * timestep)) * Qblock[j * ns + i];
                        ret[idx++] = val;
                    }
				}
            }
            for(int k = 0; k < offdn; k++) {
                if (offdi->ints[k] == i) {
                    int j = offdj->ints[k];
                    val = ( - 1 / (sigma * sigma * timestep * timestep)) * QfT[j * ns + i];
                    if (!isfinite(QfT[j * ns + i]) & debug > 0) {
                        printf("NaN in QfT\n");
                    }
                    ret[idx++] = val;
                }
			}
		}

        //blocks 1 to nt-1
        for (int t = 1; t < nt - 1; t++) {
			//calc f(T+1) = CinvG + diag(a) - 1/timestep*I
            double* fTplus1 = calloc(ns * ns, sizeof(double));
            for (int k = 0; k < CinvG->n; k++) {
                int i = CinvG->i[k];
                int j = CinvG->j[k];
                double v = CinvG->x[k];

                fTplus1[j * ns + i] = move_const * v;
                
            }
            for (int i = 0; i < ns; i++) {
                fTplus1[i * ns + i] += a_array[(t+1) * ns + i] + 1.0 / timestep;
            }
			
			
            //calc Qblock*f(T+1) and store in QfTplus1
            double* QfTplus1 = calloc(ns * ns, sizeof(double));
            if (debug > 0) printf("dgemm step");
            dgemm_("N", "N",
                &ns, &ns, &ns,
                &one,
                Qblock, &ns,
                fTplus1, &ns,
                &zero,
                QfTplus1, &ns);
            
            //calc trans(fT)*QfT and store in fTQfT
			double* fTQfT = calloc(ns * ns, sizeof(double));
            char transfT = 'T';
            if (debug > 0) printf("dgemm step");
            dgemm_(&transfT, &transB, &ns, &ns, &ns, &one, fT, &ns, QfT, &ns, &zero, fTQfT, &ns);
                
           
			

            //fill in ret for block t in order of GRAPH
             
            for (int i = t * ns; i < (t + 1) * ns; i++) {
                for (int k = 0; k < diagn; k++) {
                    if (diagi->ints[k] == i - t * ns) {
                        int j = diagj->ints[k];
                        if (j >= i - t * ns) {
							ret[idx++] = (1 / (sigma * sigma * timestep)) * fTQfT[j * ns + i - t * ns] + (1 / (sigma * sigma * timestep * timestep * timestep)) * Qblock[j * ns + i - t * ns];
                        }
                    }
                }
                for (int k = 0; k < offdn; k++) {
                    if (offdi->ints[k] == i - t * ns) {
                        int j = offdj->ints[k];
                        ret[idx++] = ( - 1 / (sigma * sigma*timestep * timestep)) * QfTplus1[j * ns + i - t * ns];
                        if (!isfinite(QfTplus1[j * ns + i]) & debug > 0) {
                            printf("NaN in QfT\n");
                        }
                    }
				}
            }
            
			//copy f(T+1) to fT for next iteration
			double* temp = fT;
			fT = fTplus1;
			fTplus1 = temp;
            double* temp2 = QfT;
            QfT = QfTplus1;
            QfTplus1 = temp2;
            free(QfTplus1);
			free(fTplus1);
        }
		//final block nt-1
        //calc trans(fT)*QfT and store in fTQfT
        double* fTQfT = calloc(ns * ns, sizeof(double));
        char transfT = 'T';
        if (debug > 0) printf("dgemm step");
        dgemm_(&transfT, &transB, &ns, &ns, &ns, &one, fT, &ns, QfT, &ns, &zero, fTQfT, &ns);

        for(int i = (nt-1)*ns; i < nt*ns; i++) {
            for (int k = 0; k < diagn; k++) {
                if (diagi->ints[k] == i - (nt - 1) * ns) {
                    int j = diagj->ints[k];
                    if (j >= i - (nt - 1) * ns) {
                        ret[idx++] = (1/(sigma * sigma * timestep)) * fTQfT[j * ns + i - (nt - 1) * ns] + (1 / (sigma * sigma * timestep * timestep * timestep)) * Qblock[j * ns + i - (nt - 1) * ns];
                    }
                }
            }
		}

        /*for (int i = (nt - 1) * ns; i < nt * ns; i++) {
            for (int j = i; j < nt * ns; j++) {
                ret[idx++] = (sigma * sigma / timestep) * fTQfT[(j - (nt - 1) * ns) * ns + (i - (nt - 1) * ns)];
            }
		}*/
        
        free(Qblock);
        free(fT);
        free(QfT);
		free(fTQfT);
        
        assert(idx == M + 2);

    }
    break;
    case INLA_CGENERIC_MU:
    {
        // return (N, mu)
        if (debug > 0) {
            printf("INLA_CGENERIC_MU\n");
        }
        ret = calloc(1 + N, sizeof(double));
        assert(ret);
        ret[0] = N; /* dimension */

        inla_cgeneric_mat_tp* L_mat = calloc(1, sizeof(inla_cgeneric_mat_tp));
        L_mat->x = calloc(N * N, sizeof(double));
        L_mat->nrow = N;
        L_mat->ncol = N;
        Lmat(growth, carry_cap, move_const, timestep, linpoint->doubles, ns, nt, CinvG, L_mat->x);

        inla_cgeneric_vec_tp* rvector = calloc(1, sizeof(inla_cgeneric_vec_tp));
        rvector->doubles = calloc(N, sizeof(double));
        rvector->len = N;
        r_vector(growth, carry_cap, move_const, linpoint->doubles, mag_grad_sq->doubles, ns, nt, rvector->doubles);
        for (int i = 0; i < ns; i++) {
            rvector->doubles[i] = prior_mean->doubles[i];
        }

        //calculate L_mat^-1 * rvector
        int* ipiv = calloc(ns * nt, sizeof(int));
        int lda = N;
        int ldb = N;
        int nrhs = 1;
        int info;
        double* A = calloc(N * N,sizeof(double));
        double* B = calloc(N, sizeof(double));
		memcpy(A, L_mat->x, N* N * sizeof(double));
        memcpy(B, rvector->doubles, N * sizeof(double));
        if (debug > 0) printf("dgesv step");
        dgesv_(&N, &nrhs, A, &lda, ipiv, B, &ldb, &info);
        if (info != 0) {
            printf("dgesv failed, info = %d\n", info);
            free(A); free(B); free(ipiv);
            free(L_mat->x); free(L_mat);
            free(rvector->doubles); free(rvector);
            return NULL;
        }
        for (int i = 0; i < N; i++) {
            if (!isfinite(B[i])) {
                printf("NaN in MU solution at %d\n", i);
            }
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
            normal_pdf_log(theta[2], pmove->doubles[0], pmove->doubles[1]) +
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
