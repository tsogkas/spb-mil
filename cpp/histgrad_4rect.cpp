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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 4)
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs != 3)
        mexErrMsgTxt("Wrong number of output arguments");

    /* Input arguments */
    unsigned short int *im = (unsigned short int*) mxGetData(prhs[0]);
    const int nbins = mxGetScalar(prhs[1]);
    const int scale = mxGetScalar(prhs[2]);
    const int ratio = mxGetScalar(prhs[3]);
    const int nrows = mxGetM(prhs[0]);
    const int ncols = mxGetN(prhs[0]);

    /* Output arguments (tg_dlc, tg_drc, tg_dlr) */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    float *diffup   = (float*) mxGetData(plhs[0]);
    float *diffdown = (float*) mxGetData(plhs[1]);
    float *diffud   = (float*) mxGetData(plhs[2]);

    /* Internal arguments for function */
    float *iir, tgC1, tgC2, tgU, tgD;
    const float EPS  = 0.00000001;
    const int halfb  = ratio*scale;


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
        const int jstart = halfb+1;
        const int jend   = ncols-halfb;
        const int istart = 2*scale+2;
        const int iend   = nrows-2*scale-1;
        const int temp1  = 2*scale+2;
        const int temp2  = 2*scale+1;
        int temp3 = 0, temp4 = 0;
        for (int j=jstart; j<jend; ++j) {
            temp3 = nrows*(j-halfb-1);
            temp4 = nrows*(j+halfb);
            for (int i=istart; i<iend; ++i) {
                // Calculate tgC, tgL, tgR
                tgU = (iir[(i-temp1)+temp3] -
                       iir[(i-temp1)+temp4] -
                       iir[(i-scale-1)+temp3] +
                       iir[(i-scale-1)+temp4]);
                tgC1= (iir[(i-scale-1)+temp3] -
                       iir[(i-scale-1)+temp4] -
                       iir[i+temp3] +
                       iir[i+temp4]);
                tgC2= (iir[(i-1)+temp3] -
                       iir[(i-1)+temp4] -
                       iir[(i+scale)+temp3] +
                       iir[(i+scale)+temp4]);
                tgD = (iir[(i+scale)+temp3] -
                       iir[(i+scale)+temp4] -
                       iir[(i+temp2)+temp3] +
                       iir[(i+temp2)+temp4]);
                // Calculate histogram differences
                diffup[i+nrows*j]   += ((tgU-tgC1)*(tgU-tgC1))/((tgU+tgC1+EPS));
                diffdown[i+nrows*j] += ((tgC2-tgD)*(tgC2-tgD))/((tgC2+tgD+EPS));
                diffud[i+nrows*j]   += ((tgU-tgD)*(tgU-tgD))/((tgU+tgD+EPS));
            }
        }
    }
    mxFree(iir);
}
