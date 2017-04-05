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
    /* Check number of input arguments */
    if (nrhs != 4)
        mexErrMsgTxt("Wrong number of input arguments");

    /* Check number of output arguments */
    if (nlhs != 3)
        mexErrMsgTxt("Wrong number of output arguments");

    /* Input arguments */
    const unsigned short int *im = (const unsigned short int*) mxGetData(prhs[0]);
    const int nbins = (const int) mxGetScalar(prhs[1]);
    const int scale = (const int) mxGetScalar(prhs[2]);
    const int ratio = (const int) mxGetScalar(prhs[3]);
    const int nrows = (const int) mxGetM(prhs[0]);
    const int ncols = (const int) mxGetN(prhs[0]);

    /* Output arguments (tg_dlc, tg_drc, tg_dlr) */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    float *dlc = (float*) mxGetData(plhs[0]);
    float *drc = (float*) mxGetData(plhs[1]);
    float *dlr = (float*) mxGetData(plhs[2]);

    float *iir	     = (float*) mxMalloc(nrows*ncols*sizeof(float));
    const float EPS  = 0.00000001;
    const int halfb  = ratio*scale;
    const int istart = 3*scale+2;		// set image limits (do only necessary computations)
    const int istep  = 2*scale+1;
    const int iend   = istart+istep;
    const int jstart = halfb+1;
    const int jend   = ncols-halfb;
    const int wstop  = nrows-3*scale-1;
    const int temp1  = 3*scale+1;
    const int temp2  = 3*scale+2;

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
        int old = 0, next = 0, ci = 0, temp3 = 0, temp4 = 0;
        float otgC, otgU, otgD, ntgD;
//#pragma omp parallel for private(old,next,ci,temp3,temp4,otgL,otgC,otgR,ntgR)
        for (int j=jstart; j<jend; ++j) {
            temp3 = nrows*(j-halfb-1);	// left rectangle nodes
            temp4 = nrows*(j+halfb);    // right rectangle nodes
            for (int i=istart; i<iend; ++i)	{
                old = i+nrows*j;
                ci = i;
                otgU = (iir[(i-temp2)+temp3]   	-   // upper rectangle
                        iir[(i-temp2)+temp4]  	-
                        iir[(i-scale-1)+temp3]  +
                        iir[(i-scale-1)+temp4]);
                otgC = (iir[(i-scale-1)+temp3]  -   // center rectangle
                        iir[(i-scale-1)+temp4]  -
                        iir[(i+scale)+temp3]    +
                        iir[(i+scale)+temp4]);
                otgD = (iir[(i+scale)+temp3] 	-   // bottom rectangle
                        iir[(i+scale)+temp4] 	-
                        iir[(i+temp1)+temp3] 	+
                        iir[(i+temp1)+temp4]);
                dlc[old] += ((otgU-otgC)*(otgU-otgC)) / ((otgU+otgC+EPS));
                drc[old] += ((otgD-otgC)*(otgD-otgC)) / ((otgD+otgC+EPS));
                dlr[old] += ((otgU-otgD)*(otgU-otgD)) / ((otgU+otgD+EPS));

                // traverse image in a way that already computed values can be reused
                while (ci+istep < wstop) {
                    ci  += istep;
                    ntgD = (iir[(ci+scale)+temp3] -
                            iir[(ci+scale)+temp4] -
                            iir[(ci+temp1)+temp3] +
                            iir[(ci+temp1)+temp4]);
                    next = ci+nrows*j;
                    dlc[next]  = drc[old];
                    drc[next] += ((ntgD-otgD)*(ntgD-otgD)) / ((ntgD+otgD+EPS));
                    dlr[next] += ((otgC-ntgD)*(otgC-ntgD)) / ((otgC+ntgD+EPS));
                    otgC = otgD;
                    otgD = ntgD;
                    old  = next;
                }
            }
        }
    }
    mxFree(iir);
}
