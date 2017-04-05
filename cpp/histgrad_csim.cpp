/*
 * histgrad_csim.cpp
 *
 * Calculates histogram features for symmetry detector based on color similarity,
 * using rotated versions of the original image and integral image representation.
 *
 *  Created on: Nov 12, 2012
 *      Author: Stavros Tsogkas
 */

#include <cmath>
#include "matrix.h"
#include "mex.h"


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
	plhs[0] = mxCreateNumericMatrix(nrows*ncols,nbins,mxSINGLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericMatrix(nrows*ncols,nbins,mxSINGLE_CLASS,mxREAL);
	plhs[2] = mxCreateNumericMatrix(nrows*ncols,nbins,mxSINGLE_CLASS,mxREAL);
	float *tg_dlc = (float*) mxGetData(plhs[0]);
	float *tg_drc = (float*) mxGetData(plhs[1]);
	float *tg_dlr = (float*) mxGetData(plhs[2]);

	/* Internal arguments for function */
	float *iir, otgC, otgL, otgR, ntgR;
	iir   = (float*) mxMalloc(nrows*ncols*sizeof(float));

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

		// define limits for computation
		int start = 3*scale+2;
		int istep = 2*scale+1;
		int iend = start+istep-1;	// part of the image where features are computed without reuse
		int jend = ncols-3*scale-1;
		int old,next, ci;
		// calculate histogram differences
		for (int j=start; j<=jend; ++j)
			for (int i=start; i<=iend; ++i)
			{
				old = i+nrows*j+(nrows*ncols)*(ibin-1);
				ci = i;
				// Calculate tgC, tgL, tgR
				otgL = (iir[(i-3*scale-2)+nrows*(j-3*scale-1)] - iir[(i-3*scale-2)+nrows*(j+3*scale)] - iir[(i-scale-1)+nrows*(j-3*scale-1)] + iir[(i-scale-1)+nrows*(j+3*scale)]);
				otgC = (iir[(i-scale-1)+nrows*(j-3*scale-1)] - iir[(i-scale-1)+nrows*(j+3*scale)] - iir[(i+scale)+nrows*(j-3*scale-1)] + iir[(i+scale)+nrows*(j+3*scale)]);
				otgR = (iir[(i+scale)+nrows*(j-3*scale-1)] - iir[(i+scale)+nrows*(j+3*scale)] - iir[(i+3*scale+1)+nrows*(j-3*scale-1)] + iir[(i+3*scale+1)+nrows*(j+3*scale)]);
				// Calculate histogram differences
				tg_dlc[old] = abs(otgL-otgC);
				tg_drc[old] = abs(otgR-otgC);
				tg_dlr[old] = abs(otgL-otgR);
				// traverse image in a way that already computed values can be reused
				while (ci+istep <= nrows-3*scale-1)
				{
					ci += istep;
					ntgR = (iir[(ci+scale)+nrows*(j-3*scale-1)] - iir[(ci+scale)+nrows*(j+3*scale)] - iir[(ci+3*scale+1)+nrows*(j-3*scale-1)] + iir[(ci+3*scale+1)+nrows*(j+3*scale)]);
					next = ci+nrows*j+(nrows*ncols)*(ibin-1);
					tg_dlc[next] = tg_drc[old];
					tg_drc[next] = abs(ntgR-otgR);
					tg_dlr[next] = abs(otgC-ntgR);
					otgC = otgR;
					otgR = ntgR;
					old = next;
				}

			}
	}
	mxFree(iir);
}


