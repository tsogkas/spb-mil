/*
 * histgrad2.cpp
 *
 * Calculates the histogram differences for symmetry detector using rotated versions of
 * the original image and integral image representation.
 *
 *  Latest revision in: April 2013
 *      Author: Stavros Tsogkas
 *		<stavros.tsogkas[at]ecp.fr>
 */

#include <cmath>
#include "matrix.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 4)
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of output arguments");

    /* Input arguments */
    const unsigned short int *im = (const unsigned short int*) mxGetData(prhs[0]);
    const int nbins = (const int) mxGetScalar(prhs[1]);
    const int scale = (const int) mxGetScalar(prhs[2]);
    const int ratio = (const int) mxGetScalar(prhs[3]);
    const int nrows = (const int) mxGetM(prhs[0]);
    const int ncols = (const int) mxGetN(prhs[0]);

    /* Output arguments  */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    float *dlc = (float*) mxGetData(plhs[0]);

    float otgC, otgL;
    float *iir	    = (float*) mxMalloc(nrows*ncols*sizeof(float));
    const float EPS = 0.00000001;
    const int halfb = ratio*scale;
    const int start = halfb+2;		 // image limits (do only necessary computations)
    const int iend  = nrows-halfb-1; // part of the image where features are computed without reuse
    const int jend  = ncols-halfb-1;
    const int temp1 = halfb+1;
    const int temp2 = halfb-2;


    for (int ibin=1; ibin<=nbins; ++ibin) {
        // Integral image representation
        for (int i=0; i<nrows; ++i)
            iir[i] = (im[i]==ibin ? 1 : 0);
        for (int j=1; j<ncols; ++j) 	// cumulative sum over rows
            for (int i=0; i<nrows; ++i)
                iir[i+nrows*j] = iir[i+nrows*(j-1)] + (im[i+nrows*j]==ibin ? 1 : 0);
        for (int j=0; j<ncols; ++j) 	// cumulative sum over columns
            for (int i=1; i<nrows; ++i)
                iir[i+nrows*j] += iir[i+nrows*j-1];

        // Calculate histogram differences
        int temp3 = 0, temp4 = 0;
        for (int j=start; j<=jend; ++j) {
            temp3 = nrows*(j-halfb-1);	// cache calculated indexes for speed
            temp4 = nrows*(j+halfb);
            for (int i=start; i<=iend; ++i)	{
                // Calculate tgC, tgL
                otgL = (iir[(i-temp2)+temp3]   	-
                        iir[(i-temp2)+temp4]  	-
                        iir[(i-scale-1)+temp3]  +
                        iir[(i-scale-1)+temp4]);
                otgC = (iir[(i-scale-1)+temp3]  -
                        iir[(i-scale-1)+temp4]  -
                        iir[(i+scale)+temp3]    +
                        iir[(i+scale)+temp4]);
                dlc[i+nrows*j] += ((otgL-otgC)*(otgL-otgC)) / ((otgL+otgC+EPS));
            }
        }
    }
    mxFree(iir);
}
