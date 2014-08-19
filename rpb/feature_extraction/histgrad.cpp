/*
 * histgrad.cpp
 *
 * Calculates the histogram differences for symmetry detector using rotated versions of
 * the original image and integral image representation.
 *
 *  Created on: Apr 11, 2012
 *      Author: Stavros Tsogkas
 */

#include <cmath>
#include "matrix.h"
#include "mex.h"
#include <omp.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Check number of input arguments */
	if (nrhs != 3)
		mexErrMsgTxt("Three input arguments are required (input image, number of bins and scale).");

	/* Check number of output arguments */
	if (nlhs != 3)
		mexErrMsgTxt("Three output arguments are required (the three histogram differences).");

	/* Input arguments */
	unsigned short int *im = (unsigned short int*) mxGetData(prhs[0]);
	int nbins = mxGetScalar(prhs[1]);
	int scale = mxGetScalar(prhs[2]);
	int nrows = mxGetM(prhs[0]);
	int ncols = mxGetN(prhs[0]);

	/* Output arguments (tg_dlc, tg_drc, tg_dlr) */
	plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
	plhs[2] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
	float *tg_dlc = (float*) mxGetData(plhs[0]);
	float *tg_drc = (float*) mxGetData(plhs[1]);
	float *tg_dlr = (float*) mxGetData(plhs[2]);

	/* Internal arguments for function */
	float *iir, tgC, tgL, tgR;
	const float EPS = 0.00000001, rectArea = (2*scale+1)*(6*scale+1);
	iir   = (float*) mxMalloc(nrows*ncols*sizeof(float));

	//#pragma omp parallel for
	for (int ibin=1; ibin<=nbins; ++ibin)
	{
		// integral image representation
		for (int i=0; i<nrows; ++i)
			iir[i] = (im[i]==ibin ? 1 : 0);
		for (int j=1; j<ncols; ++j)	// cumulative sum over rows
			for (int i=0; i<nrows; ++i)
				iir[i+nrows*j] = iir[i+nrows*(j-1)] + (im[i+nrows*j]==ibin ? 1 : 0);
		for (int j=0; j<ncols; ++j)	// cumulative sum over columns
			for (int i=1; i<nrows; ++i)
				iir[i+nrows*j] += iir[i-1+nrows*j];

		// calculate histogram differences
		for (int j=3*scale+2, jend=ncols-3*scale-1; j<=jend; ++j)
			for (int i=3*scale+2, iend=nrows-3*scale-1; i<=iend; ++i)
			{
				// Calculate tgC, tgL, tgR
				tgL = (iir[(i-3*scale-2)+nrows*(j-3*scale-1)] - iir[(i-3*scale-2)+nrows*(j+3*scale)] - iir[(i-scale-1)+nrows*(j-3*scale-1)] + iir[(i-scale-1)+nrows*(j+3*scale)]) / rectArea;
				tgC = (iir[(i-scale-1)+nrows*(j-3*scale-1)] - iir[(i-scale-1)+nrows*(j+3*scale)] - iir[(i+scale)+nrows*(j-3*scale-1)] + iir[(i+scale)+nrows*(j+3*scale)]) / rectArea;
				tgR = (iir[(i+scale)+nrows*(j-3*scale-1)] - iir[(i+scale)+nrows*(j+3*scale)] - iir[(i+3*scale+1)+nrows*(j-3*scale-1)] + iir[(i+3*scale+1)+nrows*(j+3*scale)]) / rectArea;
				// Calculate histogram differences
				tg_dlc[i+nrows*j] += ((tgL-tgC)*(tgL-tgC))/((tgL+tgC+EPS));
				tg_drc[i+nrows*j] += ((tgR-tgC)*(tgR-tgC))/((tgR+tgC+EPS));
				tg_dlr[i+nrows*j] += ((tgL-tgR)*(tgL-tgR))/((tgL+tgR+EPS));
			}
	}
	mxFree(iir);
}
