//***************************************************************************
//
// Matlab C routine file:  anaskel.cpp
//
// Written 8/04 by N. Howe
//
// Input:
//   skel:  shape to find skeleton for
//
// Output:
//   dmap:  directional maps
//   exy:   coordinates of endpoints
//   jxy:   coordinates of junction points
//
//***************************************************************************

#include "mex.h"

//***************************************************************************

#define SQR(x) (x)*(x)
#define MIN(x,y) (((x) < (y)) ? (x):(y))
#define MAX(x,y) (((x) > (y)) ? (x):(y))
#define ABS(x) (((x) < 0) ? (-(x)):(x))
#define errCheck(a,b) if (!(a)) mexErrMsgTxt((b));
#define bdCheck(a,b) mxAssert(((a)>=0)&&((a)<(b)),"Bounds error.")

//***************************************************************************

const int N_DIRECTIONS = 4;

const int connected_nbrs[256] = {
  0,1,1,1,1,1,1,1,
  1,2,2,2,1,1,1,1,
  1,2,2,2,1,1,1,1,
  1,2,2,2,1,1,1,1,
  1,2,2,2,2,2,2,2,
  2,3,3,3,2,2,2,2,
  1,2,2,2,1,1,1,1,
  1,2,2,2,1,1,1,1,
  1,1,2,1,2,1,2,1,
  2,2,3,2,2,1,2,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,2,1,2,1,
  2,2,3,2,2,1,2,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,2,1,2,1,
  2,2,3,2,2,1,2,1,
  2,2,3,2,2,1,2,1,
  2,2,3,2,2,1,2,1,
  2,2,3,2,3,2,3,2,
  3,3,4,3,3,2,3,2,
  2,2,3,2,2,1,2,1,
  2,2,3,2,2,1,2,1,
  1,1,2,1,2,1,2,1,
  2,2,3,2,2,1,2,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,2,1,2,1,
  2,2,3,2,2,1,2,1,
  1,1,2,1,1,1,1,1,
  1,1,2,1,1,1,1,1
};

const int nbr_branches[256] = {
  0,1,1,1,1,2,1,2,
  1,2,2,2,1,2,2,2,
  1,2,2,2,2,3,2,3,
  1,2,2,2,2,3,2,3,
  1,2,2,2,2,3,2,3,
  2,3,3,3,2,3,3,3,
  1,2,2,2,2,3,2,3,
  2,3,3,3,2,3,3,3,
  1,2,2,2,2,3,2,3,
  2,3,3,3,2,3,3,3,
  2,3,3,3,3,4,3,4,
  2,3,3,3,3,4,3,4,
  1,2,2,2,2,3,2,3,
  2,3,3,3,2,3,3,3,
  2,3,3,3,3,4,3,4,
  2,3,3,3,3,4,3,4,
  1,1,2,2,2,2,2,2,
  2,2,3,3,2,2,3,3,
  2,2,3,3,3,3,3,3,
  2,2,3,3,3,3,3,3,
  2,2,3,3,3,3,3,3,
  3,3,4,4,3,3,4,4,
  2,2,3,3,3,3,3,3,
  3,3,4,4,3,3,4,4,
  1,2,2,2,2,3,2,3,
  2,3,3,3,2,3,3,3,
  2,3,3,3,3,4,3,4,
  2,3,3,3,3,4,3,4,
  2,2,3,3,3,3,3,3,
  3,3,4,4,3,3,4,4,
  2,3,3,3,3,4,3,4,
  3,3,4,4,3,4,4,4,
};

const int isN[256] = {
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,1,1,1,0,1,1,1,
  0,1,1,1,0,1,1,1,
  0,1,0,1,0,1,0,1,
  0,1,0,1,0,1,0,1,
  0,1,1,1,0,1,1,1,
  0,1,1,1,0,1,1,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,1,1,1,0,1,1,1,
  0,1,1,1,0,1,1,1,
  0,1,0,1,0,1,0,1,
  0,1,0,1,0,1,0,1,
  0,1,1,1,0,1,1,1,
  0,1,1,1,0,1,1,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,1,0,1,0,1,0,1,
  0,1,0,1,0,1,0,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,1,0,1,0,1,0,1,
  0,1,0,1,0,1,0,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1
};

const int isNE[256] = {
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,0,1,1,0,0,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1,
  0,1,1,1,1,1,1,1
};

const int isE[256] = {
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,0,0,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1
};

const int isSE[256] = {
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,0,0,0,0,0,0,0,
  0,1,0,1,0,1,0,1,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,
  0,0,0,0,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,0,0,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,0,0,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  0,0,0,0,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1
};

//****************************************************************************
//
// neighborhood() returns a byte value representing the immediate neighborhood
// of the specified point in bit form, clockwise from 12:00.
//
// R/O arr:  binary pixel array
// R/O i,j:  coordinates of point
// R/O nrow, ncol:  dimensions of binary image
// Returns:  bit representation of 8-neighbors
//

template <class T>
inline int neighborhood(T *arr, int i, int j, int nrow, int ncol) {
  int p = i+j*nrow;
  int condition = 8*(i <= 0)+4*(j <= 0)+2*(i >= nrow-1)+(j >= ncol-1);

  //mexPrintf("Condition:  %d\n",condition);
  switch (condition) {
  case 0:  // all sides valid
    mxAssert((p-nrow-1)>=0,"Too low.");
    mxAssert((p+nrow+1)<nrow*ncol,"Too high.");
    return (arr[p-1]?1:0) + (arr[p+nrow-1]?2:0) + (arr[p+nrow]?4:0) +
      (arr[p+nrow+1]?8:0) + (arr[p+1]?16:0) + (arr[p-nrow+1]?32:0) +
      (arr[p-nrow]?64:0) + (arr[p-nrow-1]?128:0);
  case 1:  // right side not valid
    return (arr[p-1]?1:0) + (arr[p+1]?16:0) + (arr[p-nrow+1]?32:0) +
      (arr[p-nrow]?64:0) + (arr[p-nrow-1]?128:0);
  case 2:  // bottom not valid
    return (arr[p-1]?1:0) + (arr[p+nrow-1]?2:0) + (arr[p+nrow]?4:0) +
      (arr[p-nrow]?64:0) + (arr[p-nrow-1]?128:0);
  case 3:  // bottom and right not valid
    return (arr[p-1]?1:0) + (arr[p-nrow]?64:0) + (arr[p-nrow-1]?128:0);
  case 4:  // left side not valid
    return (arr[p-1]?1:0) + (arr[p+nrow-1]?2:0) + (arr[p+nrow]?4:0) +
      (arr[p+nrow+1]?8:0) + (arr[p+1]?16:0);
  case 5:  // left and right sides not valid
    return (arr[p-1]?1:0) + (arr[p+1]?16:0);
  case 6:  // left and bottom sides not valid
    return (arr[p-1]?1:0) + (arr[p+nrow-1]?2:0) + (arr[p+nrow]?4:0);
  case 7:  // left, bottom, and right sides not valid
    return (arr[p-1]?1:0);
  case 8:  // top side not valid
    return (arr[p+nrow]?4:0) + (arr[p+nrow+1]?8:0) + (arr[p+1]?16:0) + 
      (arr[p-nrow+1]?32:0) + (arr[p-nrow]?64:0);
  case 9:  // top and right not valid
    return (arr[p+1]?16:0) + (arr[p-nrow+1]?32:0) + (arr[p-nrow]?64:0);
  case 10:  // top and bottom not valid
    return (arr[p+nrow]?4:0) + (arr[p-nrow]?64:0);
  case 11:  // top, bottom and right not valid
    return (arr[p-nrow]?64:0);
  case 12:  // top and left not valid
    return (arr[p+nrow]?4:0) + (arr[p+nrow+1]?8:0) + (arr[p+1]?16:0);
  case 13:  // top, left and right sides not valid
    return (arr[p+1]?16:0);
  case 14:  // top, left and bottom sides not valid
    return (arr[p+nrow]?4:0);
  case 15:  // no sides valid
    return 0;
  }
}

//****************************************************************************
//
// dotrim() trims unnecessary points from the skeleton.
//
// R/O inp:  binary pixel array
// R/O nrow, ncol:  dimensions of binary image
// W/O skel:  newly trimmed skeleton
//

template <class T>
void dotrim(T *inp, int nrow, int ncol, unsigned char *skel) {
  int i, j, hood;

  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
	skel[i+j*nrow] = (unsigned char)(inp[i+j*nrow]);
    }
  }
  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
      if (skel[i+j*nrow]) {
	hood = neighborhood(skel,i,j,nrow,ncol);
	skel[i+j*nrow] = (connected_nbrs[hood] > 1)||(nbr_branches[hood]==1);
      } else {
	skel[i+j*nrow] = 0;
      }
    }
  }
}

//****************************************************************************
//
// This is the Matlab entry point.
//
//****************************************************************************

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i, j, p, hood;
  int nrow, ncol, nend, nbar, njunc, iend, ijunc;
  unsigned char *skel;
  double *north, *northeast, *east, *southeast, *exy, *jxy;
  mxArray *cell;
	
  // check for proper number of arguments
  errCheck(nrhs == 1,"Exactly one input argument required.");
  errCheck(nlhs <= 3,"Too many output arguments.");

  // check format of arguments
  errCheck(mxIsUint8(prhs[0])||mxIsLogical(prhs[0])||mxIsDouble(prhs[0]),
           "Input must be binary image.");
  nrow = mxGetM(prhs[0]);
  ncol = mxGetN(prhs[0]);

  // allocate temporary space
  skel = (unsigned char *)mxMalloc(nrow*ncol*sizeof(unsigned char));

  // trim skeleton
  if (mxIsDouble(prhs[0])) {
    dotrim(mxGetPr(prhs[0]),nrow,ncol,skel);
  } else {
    dotrim((unsigned char *)(mxGetPr(prhs[0])),nrow,ncol,skel);
  }  
  
  // allocate output space
  plhs[0] = mxCreateCellMatrix(1,4);
  cell = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  mxSetCell(plhs[0],0,cell);
  north = mxGetPr(cell);
  cell = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  mxSetCell(plhs[0],1,cell);
  northeast = mxGetPr(cell);
  cell = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  mxSetCell(plhs[0],2,cell);
  east = mxGetPr(cell);
  cell = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  mxSetCell(plhs[0],3,cell);
  southeast = mxGetPr(cell);

  // analyze...
  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) {
      p = i+j*nrow;
      if (skel[p]) {
	hood = neighborhood(skel,i,j,nrow,ncol);
	north[p] = isN[hood];
	northeast[p] = isNE[hood];
	east[p] = isE[hood];
	southeast[p] = isSE[hood];
	//mexPrintf("Point (%d,%d):  Hood %d -> %d, %d, %d, %d.\n",
	//	  i,j,hood,isN[hood],isNE[hood],isE[hood],isSE[hood]);
      } else {
	north[p] = northeast[p] = east[p] = southeast[p] = 0;
      }
    }
  }

  // extra data if necessary
  if (nlhs > 1) {
    // count junctions and endpoints
    njunc = 0;
    nend = 0;
    nbar = 0;
    for (j = 0; j < ncol; j++) {
      for (i = 0; i < nrow; i++) {
	if (skel[i+j*nrow]) {
	  hood = neighborhood(skel,i,j,nrow,ncol);
	  switch (nbr_branches[hood]) {
	  case 0:
	  case 1:
	    nend++;
	    break;
	  case 2:
	    nbar++;
	    break;
	  case 3:
	  case 4:
	    njunc++;
	    break;
	  }
	}
      }
    }
    //mexPrintf("Counted.\n");

    plhs[1] = mxCreateDoubleMatrix(2, nend, mxREAL);
    exy = mxGetPr(plhs[1]);
    iend = 0;
    for (j = 0; j < ncol; j++) {
      for (i = 0; i < nrow; i++) {
	if (skel[i+j*nrow]) {
	  hood = neighborhood(skel,i,j,nrow,ncol);
	  if (nbr_branches[hood] < 2) {
	    exy[iend+1] = i+1;
	    exy[iend] = j+1;
	    iend += 2;
	  }
	}
      }
    }
  }
  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(2, njunc, mxREAL);
    jxy = mxGetPr(plhs[2]);
    ijunc = 0;
    for (j = 0; j < ncol; j++) {
      for (i = 0; i < nrow; i++) {
	if (skel[i+j*nrow]) {
	  hood = neighborhood(skel,i,j,nrow,ncol);
	  if (nbr_branches[hood] > 2) {
	    jxy[ijunc+1] = i+1;
	    jxy[ijunc] = j+1;
	    ijunc += 2;
	  }
	}
      }
    }  
  }
  
  // free stuff
  mxFree(skel);
}

/****************************************************************************/
