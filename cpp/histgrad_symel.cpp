/*
 * histgrad_symel.cpp
 *
 * Calculates simple symmetry gradients (symels) using rotated versions of
 * the original image and integral image representation.
 *
 *  Latest revision in: June 2013
 *      Author: Stavros Tsogkas
 *		<stavros.tsogkas[at]ecp.fr>
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
    const unsigned short int *im = (const unsigned short int*) mxGetData(prhs[0]);
    const int nbins = (const int) mxGetScalar(prhs[1]);
    const int scale = (const int) mxGetScalar(prhs[2]);
    const int nrows = (const int) mxGetM(prhs[0]);
    const int ncols = (const int) mxGetN(prhs[0]);

    /* Output arguments (tg_dlc, tg_drc, tg_dlr) */
    plhs[0] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericMatrix(nrows,ncols,mxSINGLE_CLASS,mxREAL);
    float *tg_dlc = (float*) mxGetData(plhs[0]);
    float *tg_drc = (float*) mxGetData(plhs[1]);
    float *tg_dlr = (float*) mxGetData(plhs[2]);

    float otgC, otgL, otgR, ntgR;
    float *iir	    = (float*) mxMalloc(nrows*ncols*sizeof(float));
    const float EPS = 0.00000001;
    const int start = 3*scale+2;		// image limits (do only necessary computations)
    const int istep = 2*scale+1;
    const int iend  = start+istep-1;	// part of the image where features are computed without reuse
    const int jend  = ncols-3*scale-1;
    const int temp1 = 3*scale+1;
    const int temp2 = 3*scale-2;
    const int wstop = nrows-3*scale-1;


    for (int ibin=1; ibin<=nbins; ++ibin) {
        // Integral image representation
        for (int i=0; i<nrows; ++i)
            iir[i] = (im[i]==ibin ? 1 : 0);
        for (int j=1; j<ncols; ++j) {	// cumulative sum over rows
            const int ind1 = nrows*j;
            const int ind2 = nrows*(j-1);
            for (int i=0; i<nrows; ++i)
                iir[i+ind1] = iir[i+ind2] + (im[i+ind1]==ibin ? 1 : 0);
        }
        for (int j=0; j<ncols; ++j) {	// cumulative sum over columns
            const int ind1 = nrows*j;
            const int ind2 = nrows*j-1;
            for (int i=1; i<nrows; ++i)
                iir[i+ind1] += iir[i+ind2];
        }

        // Calculate histogram differences
        int old = 0, next = 0, ci = 0;
        for (int j=start; j<=jend; ++j) {
            const int temp3 = nrows*(j-3*scale-1);	// cache calculated indexes for speed
            const int temp4 = nrows*(j+3*scale);
            for (int i=start; i<=iend; ++i)	{
                old = i+nrows*j;
                ci = i;
                // Calculate tgC, tgL, tgR
                otgL = (iir[(i-temp2)+temp3]   	-
                        iir[(i-temp2)+temp4]  	-
                        iir[(i-scale-1)+temp3]  +
                        iir[(i-scale-1)+temp4]);
                otgC = (iir[(i-scale-1)+temp3]  -
                        iir[(i-scale-1)+temp4]  -
                        iir[(i+scale)+temp3]    +
                        iir[(i+scale)+temp4]);
                otgR = (iir[(i+scale)+temp3] 	-
                        iir[(i+scale)+temp4] 	-
                        iir[(i+temp1)+temp3] 	+
                        iir[(i+temp1)+temp4]);
                tg_dlc[old] += ((otgL-otgC)*(otgL-otgC)) / ((otgL+otgC+EPS));
                tg_drc[old] += ((otgR-otgC)*(otgR-otgC)) / ((otgR+otgC+EPS));
                tg_dlr[old] += ((otgL-otgR)*(otgL-otgR)) / ((otgL+otgR+EPS));

                // traverse image in a way that already computed values can be reused
                while (ci+istep <= wstop) {
                    ci  += istep;
                    ntgR = (iir[(ci+scale)+temp3] -
                            iir[(ci+scale)+temp4] -
                            iir[(ci+temp1)+temp3] +
                            iir[(ci+temp1)+temp4]);
                    next = ci+nrows*j;
                    tg_dlc[next]  = tg_drc[old];
                    tg_drc[next] += ((ntgR-otgR)*(ntgR-otgR)) / ((ntgR+otgR+EPS));
                    tg_dlr[next] += ((otgC-ntgR)*(otgC-ntgR)) / ((otgC+ntgR+EPS));
                    otgC = otgR;
                    otgR = ntgR;
                    old  = next;
                }

            }
        }
    }
    mxFree(iir);
}
