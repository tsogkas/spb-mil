#include "mex.h"
#include "mat.h"
#include "engine.h"
#include "math.h"

// convolution with an iir filter approximating the Gaussian
// as described in R. Deriche's '93 technical report. 
//
// Iasonas Kokkinos <jkokkin@stat.ucla.edu>

static int initialized=0;
static mxArray* intermediate_data_pos;
static mxArray* intermediate_data_neg;
static mxArray* intermediate_norm_neg;
static mxArray* intermediate_norm_pos;
static mxArray* horizontal_smoothed;

void cleanup(void)
{
	mxDestroyArray(intermediate_norm_pos);
	mxDestroyArray(intermediate_norm_neg);
	mxDestroyArray(intermediate_data_pos);
	mxDestroyArray(intermediate_data_neg);
	mxDestroyArray(horizontal_smoothed);
}

void mexFunction(int nout, mxArray *output[],int nin,  const mxArray *input[])
{
	{
		int s_m = mxGetM(input[0]);
	 	int s_n = mxGetN(input[0]);
		int i,j;
		double *pr_int_d, *pr_int_n, *pr_int_dp, *pr_int_np, *pr_int_dn, *pr_int_nn;
		double *pr_in = mxGetPr(input[0]);
		double sigma = (*mxGetPr(input[1]));
		double a_0,a_1,w_0,w_1,b_0,b_1,c_0,c_1;
        double np[5],dp[5],nn[5],dn[5];
		double accumulator, accumulator_norm;
		double sum_incoming, sum_incoming_norm;
		double sum_incoming_default;	
		int inner,offset;
        double *pr_out, *pr_dx, *pr_dxx, *pr_dxy, *pr_dyy, *pr_dy;
        
		output[0] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		output[1] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		output[2] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		output[3] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		output[4] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		output[5] = mxCreateDoubleMatrix(s_m,s_n,mxREAL);

        pr_out = mxGetPr(output[0]);
	pr_dx  = mxGetPr(output[1]);
	pr_dxx = mxGetPr(output[2]);
	pr_dy  = mxGetPr(output[3]);
	pr_dyy = mxGetPr(output[4]);
        pr_dxy = mxGetPr(output[5]);

        
		mexAtExit(cleanup);
		if (initialized==1)
		{
			int s_m_old = mxGetM(intermediate_norm_pos);
			int s_n_old = mxGetN(intermediate_norm_pos);
			if ((s_m_old!=s_m)||(s_n_old !=s_n ))
			{
				cleanup();
				initialized =0;
			}
		}
	
		if (!(initialized==1))
		{
		 intermediate_norm_pos = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
  		 intermediate_data_pos = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		 intermediate_data_neg = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		 intermediate_norm_neg = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
		 horizontal_smoothed   = mxCreateDoubleMatrix(s_m,s_n,mxREAL);
			
		 mexMakeArrayPersistent(intermediate_data_pos);
		 mexMakeArrayPersistent(intermediate_norm_pos);
		 mexMakeArrayPersistent(intermediate_data_neg);
		 mexMakeArrayPersistent(intermediate_norm_neg);
		 mexMakeArrayPersistent(horizontal_smoothed);
		 initialized = 1;
		}
        
		// SETUP DERICHE FILTERS
		a_0= 1.68;
        a_1 = 3.735;
 		w_0 = .6318;
        w_1 = 1.997;
        b_0 = 1.783;
		b_1 = 1.723;
        c_0 = -.6803;
        c_1 = -.2598;
		
		np[0] = a_0 + c_0;
		np[1] = exp(-b_1/sigma)*(c_1*sin(w_1/sigma) - (c_0 + 2*a_0)*cos(w_1/sigma)) + \
				exp(-b_0/sigma)*(a_1*sin(w_0/sigma) - (2*c_0 + a_0)*cos(w_0/sigma));
		np[2] = 2*exp(-b_0/sigma - b_1/sigma)*((a_0+c_0)*cos(w_1/sigma)*cos(w_0/sigma) -\
				a_1*cos(w_1/sigma)*sin(w_0/sigma) - cos(w_0/sigma)*c_1*sin(w_1/sigma))\
				+ c_0*exp(-2*b_0/sigma) + a_0*exp(-2*b_1/sigma);
		np[3] = exp(-b_1/sigma - 2*b_0/sigma)*(c_1*sin(w_1/sigma) - c_0*cos(w_1/sigma)) +\
				exp(-b_0/sigma - 2*b_1/sigma)*(a_1*sin(w_0/sigma) - a_0*cos(w_0/sigma));
		np[4] = 0;

		dp[4] = exp(-2*b_0/sigma - 2*b_1/sigma);
		dp[3] = -2*cos(w_0/sigma)*exp(-b_0/sigma-2*b_1/sigma) - 2*cos(w_1/sigma)*exp(-b_1/sigma-2*b_0/sigma);
		dp[2] = 4*cos(w_1/sigma)*cos(w_0/sigma)*exp(-b_0/sigma - b_1/sigma) + exp(-2*b_1/sigma) + exp(-2*b_0/sigma);
		dp[1] = -2*exp(-b_1/sigma)*cos(w_1/sigma) - 2*exp(-b_0/sigma)*cos(w_0/sigma);
		dp[0] = 1;

	    for (i=0; i<4;i++)
		{
			dn[i] = dp[i];
			nn[i] = np[i] - dp[i]*np[0];
		}
		dn[4] = dp[4];
		nn[4] = - dp[4]*np[0];

		// APPLY HORIZONTAL FILTERS
		//locator : j*size_m + i

		pr_int_d = mxGetPr(intermediate_data_pos);
		pr_int_n = mxGetPr(intermediate_norm_pos);
		sum_incoming_default =np[0] +  np[1] + np[2] + np[3] + np[4];
		// scan each matrix line 
		for (i =0; i<s_m;i++)
		{
			// and apply the two iir filters of deriche 
			for (j=0;j<4;j++)
			{
				accumulator = 0;
				accumulator_norm = 0;
				sum_incoming = 0;
				sum_incoming_norm = 0;
				offset = j*s_m +i;

				for (inner = 1; inner<=j; inner++)
				{
					accumulator = accumulator  + dp[inner]*pr_int_d[offset - inner*s_m];
					accumulator_norm = accumulator_norm + dp[inner]*pr_int_n[offset  - inner*s_m];
				}
				for (inner = 0; inner<=j; inner++)
				{
					sum_incoming  = sum_incoming + np[inner]*pr_in[offset - inner*s_m];
					sum_incoming_norm = sum_incoming_norm + np[inner];
				}
				pr_int_d[j*s_m + i] = sum_incoming - accumulator;
				pr_int_n[j*s_m + i] = sum_incoming_norm - accumulator_norm;
			}
			
			for (j=4; j<s_n; j++)
			{
				offset = j*s_m+i;
				accumulator =  -dp[1]*pr_int_d[offset  - s_m] -dp[2]*pr_int_d[offset - 2*s_m]  \
						   	   -dp[3]*pr_int_d[offset - 3*s_m] -dp[4]*pr_int_d[offset - 4*s_m];
				accumulator_norm =  -dp[1]*pr_int_n[offset  - s_m] -dp[2]*pr_int_n[offset - 2*s_m]  \
						   			-dp[3]*pr_int_n[offset - 3*s_m] -dp[4]*pr_int_n[offset - 4*s_m];

				sum_incoming  = np[0]*pr_in[offset] + np[1]*pr_in[offset - s_m] + np[2]*pr_in[offset - 2*s_m] +\
								np[3]*pr_in[offset -3*s_m] + np[4]*pr_in[offset - 4*s_m];

				pr_int_d[offset] = sum_incoming + accumulator;
				pr_int_n[offset] = sum_incoming_default + accumulator_norm;
			}
		}
		pr_int_d = mxGetPr(intermediate_data_neg);
		pr_int_n = mxGetPr(intermediate_norm_neg);
		sum_incoming_default =nn[0] +  nn[1] + nn[2] + nn[3] + nn[4];
		// apply backward filter
		for (i =0; i<s_m;i++)
		{
			for (j=0;j<4;j++)
			{
				accumulator = 0;
				accumulator_norm = 0;
				sum_incoming = 0;
				sum_incoming_norm = 0;
				offset = (s_n - j -1)*s_m +i;

				for (inner = 1; inner<=j; inner++)
				{
					accumulator = accumulator  + dn[inner]*pr_int_d[offset + inner*s_m];
					accumulator_norm = accumulator_norm + dn[inner]*pr_int_n[offset  + inner*s_m];
				}
				for (inner = 0; inner<=j; inner++)
				{
					sum_incoming  = sum_incoming + nn[inner]*pr_in[offset + inner*s_m];
					sum_incoming_norm = sum_incoming_norm + nn[inner];
				}
				pr_int_d[offset] = sum_incoming - accumulator;
				pr_int_n[offset] = sum_incoming_norm - accumulator_norm;
			}
			
			for (j=4; j<s_n; j++)
			{
				offset = (s_n - j-1)*s_m +i;
				
				accumulator =  -dn[1]*pr_int_d[offset  + s_m] -dn[2]*pr_int_d[offset + 2*s_m]  \
						   	   -dn[3]*pr_int_d[offset + 3*s_m] -dn[4]*pr_int_d[offset + 4*s_m];
				accumulator_norm =  -dn[1]*pr_int_n[offset  + s_m] -dn[2]*pr_int_n[offset + 2*s_m]  \
						   			-dn[3]*pr_int_n[offset + 3*s_m] -dn[4]*pr_int_n[offset + 4*s_m];

				sum_incoming  =  nn[1]*pr_in[offset + s_m] + nn[2]*pr_in[offset + 2*s_m] +\
								 nn[3]*pr_in[offset + 3*s_m] + nn[4]*pr_in[offset + 4*s_m];

				pr_int_d[offset] = sum_incoming + accumulator;
				pr_int_n[offset] = sum_incoming_default + accumulator_norm;
			}
		}
		pr_int_dp = mxGetPr(intermediate_data_pos);
		pr_int_np = mxGetPr(intermediate_norm_pos);
		pr_int_dn = mxGetPr(intermediate_data_neg);
		pr_int_nn = mxGetPr(intermediate_norm_neg);
		pr_in = mxGetPr(horizontal_smoothed);
		for (i =0;i<s_m;i++)
			for (j=0;j<s_n;j++)
			{
				offset = j*s_m + i;
				pr_in[offset]  = (pr_int_dp[offset] + pr_int_dn[offset])/(pr_int_np[offset] + pr_int_nn[offset]);
			}
			

		// -------------------------------------------------------------------------------	
		// NOW APPLY VERTICALLY
		// -------------------------------------------------------------------------------
		pr_int_d = mxGetPr(intermediate_data_pos);
		pr_int_n = mxGetPr(intermediate_norm_pos);
		sum_incoming_default =np[0] +  np[1] + np[2] + np[3] + np[4];
		// scan each matrix line 
		for (j =0; j<s_n;j++)
		{
			// and apply the two iir filters of deriche 
			for (i=0;i<4;i++)
			{
				accumulator = 0;
				accumulator_norm = 0;
				sum_incoming = 0;
				sum_incoming_norm = 0;
				offset = j*s_m +i;

				for (inner = 1; inner<=i; inner++)
				{
					accumulator = accumulator  + dp[inner]*pr_int_d[offset - inner];
					accumulator_norm = accumulator_norm + dp[inner]*pr_int_n[offset  - inner];
				}
				for (inner = 0; inner<=i; inner++)
				{
					sum_incoming  = sum_incoming + np[inner]*pr_in[offset - inner];
					sum_incoming_norm = sum_incoming_norm + np[inner];
				}
				pr_int_d[offset] = sum_incoming - accumulator;
				pr_int_n[offset] = sum_incoming_norm - accumulator_norm;
			}
			
			for (i=4; i<s_m; i++)
			{
				offset = j*s_m+i;
				accumulator =  -dp[1]*pr_int_d[offset - 1] -dp[2]*pr_int_d[offset - 2]  \
						   	   -dp[3]*pr_int_d[offset - 3] -dp[4]*pr_int_d[offset - 4];
				accumulator_norm =  -dp[1]*pr_int_n[offset  -1] -dp[2]*pr_int_n[offset - 2]  \
						   			-dp[3]*pr_int_n[offset - 3] -dp[4]*pr_int_n[offset - 4];

				sum_incoming  = np[0]*pr_in[offset] + np[1]*pr_in[offset - 1] + np[2]*pr_in[offset - 2] +\
								np[3]*pr_in[offset -3] + np[4]*pr_in[offset - 4];

				pr_int_d[offset] = sum_incoming + accumulator;
				pr_int_n[offset] = sum_incoming_default + accumulator_norm;
			}
		}
		
		pr_int_d = mxGetPr(intermediate_data_neg);
 		pr_int_n = mxGetPr(intermediate_norm_neg);
		sum_incoming_default =nn[0] +  nn[1] + nn[2] + nn[3] + nn[4];
		
		// apply backward filter
		for (j =0; j<s_n;j++)
		{
			for (i=0;i<4;i++)
			{
				accumulator = 0;
				accumulator_norm = 0;
				sum_incoming = 0;
				sum_incoming_norm = 0;
				offset = (j)*s_m +(s_m - i -1);

				for (inner = 1; inner<=i; inner++)
				{
					accumulator = accumulator  + dn[inner]*pr_int_d[offset + inner];
					accumulator_norm = accumulator_norm + dn[inner]*pr_int_n[offset  + inner];
				}
				for (inner = 0; inner<=i; inner++)
				{
					sum_incoming  = sum_incoming + nn[inner]*pr_in[offset + inner];
					sum_incoming_norm = sum_incoming_norm + nn[inner];
				}
				pr_int_d[offset] = sum_incoming - accumulator;
				pr_int_n[offset] = sum_incoming_norm - accumulator_norm;
			}
			
			for (i=4; i<s_m; i++)
			{
				offset = (j)*s_m +(s_m - i -1);
				
				accumulator =  -dn[1]*pr_int_d[offset + 1] -dn[2]*pr_int_d[offset + 2]  \
						   	   -dn[3]*pr_int_d[offset + 3] -dn[4]*pr_int_d[offset + 4];
				accumulator_norm =  -dn[1]*pr_int_n[offset + 1] - dn[2]*pr_int_n[offset + 2]  \
						   			-dn[3]*pr_int_n[offset + 3] - dn[4]*pr_int_n[offset + 4];

				sum_incoming  =  nn[1]*pr_in[offset + 1] + nn[2]*pr_in[offset + 2] +\
								 nn[3]*pr_in[offset + 3] + nn[4]*pr_in[offset + 4];

				pr_int_d[offset] = sum_incoming + accumulator;
				pr_int_n[offset] = sum_incoming_default + accumulator_norm;
			}
		}


 		pr_int_dp = mxGetPr(intermediate_data_pos);
		pr_int_np = mxGetPr(intermediate_norm_pos);
		pr_int_dn = mxGetPr(intermediate_data_neg);
   		pr_int_nn = mxGetPr(intermediate_norm_neg);

        
        for (i =0;i<s_m;i++)
			for (j=0;j<s_n;j++)
			{
				offset = j*s_m + i;
				pr_out[offset]  = (pr_int_dp[offset] + pr_int_dn[offset])/(pr_int_np[offset] + pr_int_nn[offset]); // + pr_int_dn[offset]);
			}

		for (i=0;i<s_m;i++)
		{
			for (j=1;j<s_n-1;j++)
			{
				offset = j*s_m +i;
				pr_dx[offset] = (pr_out[offset+s_m] - pr_out[offset-s_m])/2;
			}
			pr_dx[i] = 0; 
			pr_dx[(s_n-1)*s_m +i] = 0;
		}

		for (i=0;i<s_m;i++)
		{
			for (j=1;j<s_n-1;j++)
			{
				offset = j*s_m +i;
				pr_dxx[offset] = pr_out[offset+s_m] + pr_out[offset-s_m] - 2*pr_out[offset];
			}
			pr_dxx[i] = 0; 
			pr_dxx[(s_n-1)*s_m +i] = 0;
		}

		for (j=0;j<s_n;j++)
		{
			for (i=1;i<s_m-1;i++)
			{
				offset = j*s_m +i;
				pr_dy[offset] = (pr_out[offset+1] - pr_out[offset-1])/2;
			}
			pr_dy[j*s_m] = 0; 
			pr_dy[j*s_m + (s_m-1)] = 0;
		}

		for (j=0;j<s_n;j++)
		{
			for (i=1;i<s_m-1;i++)
			{
				offset = j*s_m +i;
				pr_dyy[offset] = pr_out[offset+1] + pr_out[offset-1] -2*pr_out[offset];
			}
			pr_dyy[j*s_m] = 0; 
			pr_dyy[j*s_m + (s_m-1)] = 0;
		}
		// remember to change the boundary

		for (j=0;j<s_n;j++)
		{
			for (i=1;i<s_m-1;i++)
			{
				offset = j*s_m +i;
				pr_dxy[offset] = (pr_dx[offset+1] - pr_dx[offset-1])/2;
			}
			pr_dxy[j*s_m] = 0; 
			pr_dxy[j*s_m + (s_m-1)] = 0;
		}
	}
}
			   
//*/