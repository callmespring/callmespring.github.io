/* some C subroutines in simulation for testing high dimensional qualitative treatment effects: 
Inference for Gaussian process using approximated restricted likelihood */

// include libraries
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<R.h>
#include<Rinternals.h>
#include<R_ext/Rdynload.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_bspline.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_statistics.h>

extern "C"{

// Fit B-spline
// Z: the set of covariates
// Y: pseduo response
// Ytest: the predicted response
// Itrain: the index of training samples
// Itest: the index of testing samples
// nnod: number of nods in spline
// ntrain: numer of training samples 
// ntest: number of testing samples
// p: dimension of covariates
void BS(double *Ytest, double *Z, double *Y, int *Itrain, int *Itest, int *nnod, int *ntest, int *ntrain, int *p) //, double *bp)
{
	try{
		/*initialize*/
		int n = (*ntest)+(*ntrain);
		gsl_bspline_workspace *gsl_bs = gsl_bspline_alloc(4, *nnod+2);

		/*use quantiles as the breakpoints*/
		gsl_matrix *gsl_bp = gsl_matrix_alloc(*nnod+2, *p);
		gsl_matrix *gsl_Z2 = gsl_matrix_alloc(*p, n);
		gsl_matrix_view gsl_Z = gsl_matrix_view_array(Z, n, *p);
		gsl_matrix_transpose_memcpy(gsl_Z2, &gsl_Z.matrix);
		double *Zvec_ptr = gsl_matrix_ptr(gsl_Z2, 0, 0);
		gsl_vector_view gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
		for (int j=0;j<*p;j++){
			gsl_sort(&Zvec_ptr[j*n], 1, n);
		}
		
		for (int j=0;j<*p;j++){
			gsl_matrix_set(gsl_bp, 0, j, Zvec_ptr[j*n]);
			for (int i=0;i<(*nnod);i++){
				gsl_matrix_set(gsl_bp, i+1, j, gsl_stats_quantile_from_sorted_data(&Zvec_ptr[j*n], 1, n, (1.0*(i+1))/(*nnod+1)));
			}
			gsl_matrix_set(gsl_bp, *nnod+1, j, Zvec_ptr[(j+1)*n-1]);
		}
		
//		gsl_matrix_view BP = gsl_matrix_view_array(bp, n, *p);
//		gsl_matrix_transpose_memcpy(&BP.matrix, gsl_Z2);
		
		gsl_matrix_free(gsl_Z2);

		/*construct the fit matrix X*/
		int ncoef = (*nnod+4)*(*p)-(*p-1);
		gsl_vector *gsl_B = gsl_vector_alloc(*nnod+4);
		gsl_matrix *gsl_X = gsl_matrix_alloc(*ntrain, ncoef);
		gsl_matrix *gsl_Cov = gsl_matrix_alloc(ncoef, ncoef);
		
		/*construct the fit matrix X[,1:(*nnod+4)]*/
		gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
		gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
		for (int i=0;i<*ntrain;i++){
			/* compute B_j(xi) for all j */
			gsl_bspline_eval(Z[(Itrain[i])*(*p)], gsl_B, gsl_bs);
			/* fill in row i of X */
			for (int k=0;k<(*nnod+4);k++){
				gsl_matrix_set(gsl_X, i, k, gsl_vector_get(gsl_B, k));
			}
		}
		
 		/*construct the fit matrix X[,(*nnod+5):end]*/
 		for (int j=1;j<*p;j++){
			gsl_bpvec = gsl_matrix_column(gsl_bp, j);
			gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
			for (int i=0;i<*ntrain;i++){
				/* compute B_j(xi) for all j */
				gsl_bspline_eval(Z[(Itrain[i])*(*p)+j], gsl_B, gsl_bs);
				/* fill in row i of X */
				for (int k=1;k<(*nnod+4);k++){
					gsl_matrix_set(gsl_X, i, j*(*nnod+3)+k, gsl_vector_get(gsl_B, k));
				}
			}
		}  

		/*construct the response*/
		gsl_vector *gsl_y = gsl_vector_alloc(*ntrain);
		for (int i=0;i<*ntrain;i++){
			gsl_vector_set(gsl_y, i, Y[Itrain[i]]);
		}

		/*least square fit*/
		gsl_vector *reg_coef = gsl_vector_alloc(ncoef);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(*ntrain, ncoef);
		double chisq;
		gsl_multifit_linear(gsl_X, gsl_y, reg_coef, gsl_Cov, &chisq, mw);

		/*prediction*/
		double yerr;
		gsl_vector *gsl_x = gsl_vector_alloc(ncoef);
		for (int i=0;i<*ntest;i++){
			// construct gsl_x[1:(*nnod+4)]
			gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
			gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
			gsl_bspline_eval(Z[Itest[i]*(*p)], gsl_B, gsl_bs);
			for (int k=0;k<(*nnod+4);k++){
				gsl_vector_set(gsl_x, k, gsl_vector_get(gsl_B, k));
			}
			// construct gsl_x[(*nndo+5):end]
			for (int j=1;j<*p;j++){
				gsl_bpvec = gsl_matrix_column(gsl_bp, j);
				gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
				gsl_bspline_eval(Z[Itest[i]*(*p)+j], gsl_B, gsl_bs);
				for (int k=1;k<(*nnod+4);k++){
					gsl_vector_set(gsl_x, j*(*nnod+3)+k, gsl_vector_get(gsl_B, k));
				}
			}
			gsl_multifit_linear_est(gsl_x, reg_coef, gsl_Cov, &Ytest[Itest[i]], &yerr);
		}

		// free the allocation
		gsl_vector_free(gsl_B);
		gsl_vector_free(reg_coef);
		gsl_vector_free(gsl_x);
		gsl_matrix_free(gsl_bp);
		gsl_matrix_free(gsl_X);
		gsl_matrix_free(gsl_Cov);
		gsl_bspline_free(gsl_bs);
		gsl_vector_free(gsl_y);
		gsl_multifit_linear_free(mw);
	}
	catch (...) {
		::Rf_error( "c++ exception (unknown reason)" );
	}
}

// Fit B-spline version II
void BSII(double *Ytest, double *Ztest, double *Z, double *Y, int *nnod, int *n0, int *p)
{
	try{
		/*initialize*/
		gsl_bspline_workspace *gsl_bs = gsl_bspline_alloc(4, *nnod+2);

		/*use quantiles as the breakpoints*/
		gsl_matrix *gsl_bp = gsl_matrix_alloc(*nnod+2, *p);
		gsl_matrix *gsl_Z2 = gsl_matrix_alloc(*p, n0[0]+n0[1]);
		for (int i=0;i<n0[0];i++){
			for (int j=0;j<*p;j++){
				gsl_matrix_set(gsl_Z2, j, i, Z[i*(*p)+j]);
			}
		}
		for (int i=0;i<n0[1];i++){
			for (int j=0;j<*p;j++){
				gsl_matrix_set(gsl_Z2, j, i+(n0[0]), Ztest[i*(*p)+j]);
			}
		}
		
		gsl_vector_view gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
		double *Zvec_ptr = gsl_matrix_ptr(gsl_Z2, 0, 0);
		for (int j=0;j<*p;j++){
			gsl_sort(&Zvec_ptr[j*(n0[0]+n0[1])], 1, n0[0]+n0[1]);
		}
		
		for (int j=0;j<*p;j++){
			gsl_matrix_set(gsl_bp, 0, j, Zvec_ptr[j*(n0[0]+n0[1])]);
			for (int i=0;i<(*nnod);i++){
				gsl_matrix_set(gsl_bp, i+1, j, gsl_stats_quantile_from_sorted_data(&Zvec_ptr[j*(n0[0]+n0[1])], 1, n0[0]+n0[1], (1.0*(i+1))/(*nnod+1)));
			}
			gsl_matrix_set(gsl_bp, *nnod+1, j, Zvec_ptr[(j+1)*(n0[0]+n0[1])-1]);
		}
		
		gsl_matrix_free(gsl_Z2);

		/*construct the fit matrix X*/
		int ncoef = (*nnod+4)*(*p)-(*p-1);
		gsl_vector *gsl_B = gsl_vector_alloc(*nnod+4);
		gsl_matrix *gsl_X = gsl_matrix_alloc(n0[0], ncoef);
		gsl_matrix *gsl_Cov = gsl_matrix_alloc(ncoef, ncoef);
		
		/*construct the fit matrix X[,1:(*nnod+4)]*/
		gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
		gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
		for (int i=0;i<n0[0];i++){
			/* compute B_j(xi) for all j */
			gsl_bspline_eval(Z[(i)*(*p)], gsl_B, gsl_bs);
			/* fill in row i of X */
			for (int k=0;k<(*nnod+4);k++){
				gsl_matrix_set(gsl_X, i, k, gsl_vector_get(gsl_B, k));
			}
		}
		
		/*construct the fit matrix X[,(*nnod+5):end]*/
		for (int j=1;j<*p;j++){
			gsl_bpvec = gsl_matrix_column(gsl_bp, j);
			gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
			for (int i=0;i<n0[0];i++){
				/* compute B_j(xi) for all j */
				gsl_bspline_eval(Z[(i)*(*p)+j], gsl_B, gsl_bs);
				/* fill in row i of X */
				for (int k=1;k<(*nnod+4);k++){
					gsl_matrix_set(gsl_X, i, j*(*nnod+3)+k, gsl_vector_get(gsl_B, k));
				}
			}
		}

		/*construct the response*/
		gsl_vector_view gsl_y = gsl_vector_view_array(Y, n0[0]);

		/*do the fit*/
		gsl_vector *reg_coef = gsl_vector_alloc(ncoef);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(n0[0], ncoef);
		double chisq;
		gsl_multifit_linear(gsl_X, &gsl_y.vector, reg_coef, gsl_Cov, &chisq, mw);

		/*prediction*/
		double yerr;
		gsl_vector *gsl_x = gsl_vector_alloc(ncoef);
		for (int i=0;i<n0[1];i++){
			// construct gsl_x[1:(*nnod+4)]
			gsl_bpvec = gsl_matrix_column(gsl_bp, 0);
			gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
			gsl_bspline_eval(Ztest[(i)*(*p)], gsl_B, gsl_bs);
			for (int k=0;k<(*nnod+4);k++){
				gsl_vector_set(gsl_x, k, gsl_vector_get(gsl_B, k));
			}
			// construct gsl_x[(*nndo+5):end]
			for (int j=1;j<*p;j++){
				gsl_bpvec = gsl_matrix_column(gsl_bp, j);
				gsl_bspline_knots(&gsl_bpvec.vector, gsl_bs);
				gsl_bspline_eval(Ztest[(i)*(*p)+j], gsl_B, gsl_bs);
				for (int k=1;k<(*nnod+4);k++){
					gsl_vector_set(gsl_x, j*(*nnod+3)+k, gsl_vector_get(gsl_B, k));
				}
			}
			gsl_multifit_linear_est(gsl_x, reg_coef, gsl_Cov, &Ytest[i], &yerr);
		}

		/*free the allocation*/
		gsl_matrix_free(gsl_bp);
		gsl_vector_free(gsl_B);
		gsl_vector_free(reg_coef);
		gsl_vector_free(gsl_x);
		gsl_matrix_free(gsl_X);
		gsl_matrix_free(gsl_Cov);
		gsl_bspline_free(gsl_bs);
		gsl_multifit_linear_free(mw);
	}
	catch (...) {
		::Rf_error( "c++ exception (unknown reason)" );
	}
}

// Tune the number of knots for fitting B-spline
// number of knots up to 5
void CV_BS(double *Ytest, double *Ztest, double *Z, double *Y, int *Itrain, int *Itest, int *ntrain, int *ntest, int *n0, int *K, int *p, int *K0)
{
	try{
		gsl_vector *gsl_SSR2 = gsl_vector_calloc(6);
		double *SSR2_ptr = gsl_vector_ptr(gsl_SSR2, 0);
		
		/*The largest training number*/
		int m0 = ntrain[5];

		/*Polynomial fit*/
		int ncoef = 3*(*p)+1;
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(m0, ncoef);
		
		/*cunstruct the fit matrix and the response, do the fit*/
		gsl_matrix *gsl_X = gsl_matrix_alloc(m0, ncoef);
		gsl_matrix *gsl_Cov = gsl_matrix_alloc(ncoef, ncoef);
		gsl_vector *gsl_B = gsl_vector_alloc(ncoef);
		gsl_vector *gsl_y = gsl_vector_alloc(m0);
		gsl_vector *reg_coef = gsl_vector_alloc(ncoef);
		double yerr;
		double chisq;
		int cumsumtrain = 0;
		int cumsumtest = 0;
		for (int k=0;k<*K;k++){
			for (int i=0;i<ntrain[k];i++){
				/*construct the fit matrix*/
				gsl_matrix_set(gsl_X, i, 0, 1.0);
				for (int j=0;j<*p;j++){
					gsl_matrix_set(gsl_X, i, 1+3*j, Z[Itrain[cumsumtrain+i]*(*p)+j]);
					gsl_matrix_set(gsl_X, i, 2+3*j, Z[Itrain[cumsumtrain+i]*(*p)+j]*Z[Itrain[cumsumtrain+i]*(*p)+j]);
					gsl_matrix_set(gsl_X, i, 3+3*j, Z[Itrain[cumsumtrain+i]*(*p)+j]*Z[Itrain[cumsumtrain+i]*(*p)+j]*Z[Itrain[cumsumtrain+i]*(*p)+j]);
				}
				/*construct the response*/
				gsl_vector_set(gsl_y, i, Y[Itrain[cumsumtrain+i]]);
			}
			for (int i=ntrain[k];i<m0;i++){
				/*construct the fit matrix*/
				for (int j=0;j<(3*(*p)+1);j++){
					gsl_matrix_set(gsl_X, i, j, 0.0);
				}
			}

			/*do the fit*/
			gsl_multifit_linear(gsl_X, gsl_y, reg_coef, gsl_Cov, &chisq, mw);

			/*prediction*/
			for (int i=0;i<ntest[k];i++){
				gsl_vector_set(gsl_B, 0, 1.0);
				for (int j=0;j<*p;j++){
					gsl_vector_set(gsl_B, 1+3*j, Z[Itest[cumsumtest+i]*(*p)+j]);
					gsl_vector_set(gsl_B, 2+3*j, Z[Itest[cumsumtest+i]*(*p)+j]*Z[Itest[cumsumtest+i]*(*p)+j]);
					gsl_vector_set(gsl_B, 3+3*j, Z[Itest[cumsumtest+i]*(*p)+j]*Z[Itest[cumsumtest+i]*(*p)+j]*Z[Itest[cumsumtest+i]*(*p)+j]);
				}
				gsl_multifit_linear_est(gsl_B, reg_coef, gsl_Cov, &Ytest[Itest[cumsumtest+i]], &yerr);
				SSR2_ptr[0] = SSR2_ptr[0]+(Ytest[Itest[cumsumtest+i]]-Y[Itest[cumsumtest+i]])*(Ytest[Itest[cumsumtest+i]]-Y[Itest[cumsumtest+i]]);
			}
			
			/*number of samples*/
			cumsumtrain = cumsumtrain + ntrain[k];
			cumsumtest = cumsumtest + ntest[k];
		}

		/*spline fit, with number of knots 1,2,3,4,5*/
		for (int nnod=1;nnod<6;nnod++){
			cumsumtrain = 0;
			cumsumtest = 0;
			for (int k=0;k<*K;k++){
				BS(Ytest, Z, Y, &Itrain[cumsumtrain], &Itest[cumsumtest], &nnod, &ntest[k], &ntrain[k], p);//, bp);
				for (int i=0;i<ntest[k];i++){
					SSR2_ptr[nnod] = SSR2_ptr[nnod]+(Ytest[Itest[cumsumtest+i]]-Y[Itest[cumsumtest+i]])*(Ytest[Itest[cumsumtest+i]]-Y[Itest[cumsumtest+i]]);
				}
				cumsumtrain = cumsumtrain + ntrain[k];
				cumsumtest = cumsumtest + ntest[k];
			}
		}

		/*return the index with the smallest SSR2*/
		int MI = gsl_vector_min_index(gsl_SSR2);
		K0[0] = MI;
		if (MI==0){
			gsl_X = gsl_matrix_alloc(n0[0], ncoef);
			gsl_y = gsl_vector_alloc(n0[0]);
			mw = gsl_multifit_linear_alloc(n0[0], ncoef);

			for (int i=0;i<n0[0];i++){
				/*construct the fit matrix*/
				gsl_matrix_set(gsl_X, i, 0, 1.0);
				for (int j=0;j<*p;j++){
					gsl_matrix_set(gsl_X, i, 1+3*j, Z[i*(*p)+j]);
					gsl_matrix_set(gsl_X, i, 2+3*j, Z[i*(*p)+j]*Z[i*(*p)+j]);
					gsl_matrix_set(gsl_X, i, 3+3*j, Z[i*(*p)+j]*Z[i*(*p)+j]*Z[i*(*p)+j]);
				}
				/*construct the response*/
				gsl_vector_set(gsl_y, i, Y[i]);
			}

			/*do the fit*/
			gsl_multifit_linear(gsl_X, gsl_y, reg_coef, gsl_Cov, &chisq, mw);

			/*prediction*/
			for (int i=0;i<n0[1];i++){
				gsl_vector_set(gsl_B, 0, 1.0);
				for (int j=0;j<*p;j++){
					gsl_vector_set(gsl_B, 1+3*j, Ztest[i*(*p)+j]);
					gsl_vector_set(gsl_B, 2+3*j, Ztest[i*(*p)+j]*Ztest[i*(*p)+j]);
					gsl_vector_set(gsl_B, 3+3*j, Ztest[i*(*p)+j]*Ztest[i*(*p)+j]*Ztest[i*(*p)+j]);
				}
				gsl_multifit_linear_est(gsl_B, reg_coef, gsl_Cov, &Ytest[i], &yerr);
			}
		}
		else{
			BSII(Ytest, Ztest, Z, Y, &MI, n0, p);
		}

		/*free the allocation*/
		gsl_vector_free(gsl_SSR2);
		gsl_vector_free(gsl_y);
		gsl_vector_free(reg_coef);
		gsl_matrix_free(gsl_X);
		gsl_vector_free(gsl_B);
		gsl_matrix_free(gsl_Cov);
		gsl_multifit_linear_free(mw);
	}
	catch (...) {
		::Rf_error( "c++ exception (unknown reason)" );
	}
}

}