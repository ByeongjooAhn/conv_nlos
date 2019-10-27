#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Check for proper number of arguments */

	if (nrhs != 13) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
				"mexComputeBackprojection requires 13 input arguments.");
	}

    // INPUT
    double  *transient, *x, *y, *z, *l, *s, *n_voxel;
    double  tbin1, delta, z_wall;
    int     n_l, n_s, n_t;

    // OUTPUT
    double  *rho;

    // TEMP
    double p_x, p_y, p_z, l_x, l_y, s_x, s_y, r_l, r_s, r_l2, r_s2, r, m, lastz, idx_r;
    int idx_r_floor, idx_r_ceil, idx_zyx, idx_ls, idx_lst, firstNZ, lastNZ;

    transient   =   mxGetPr(prhs[0]);
    x           =   mxGetPr(prhs[1]);
    y           =   mxGetPr(prhs[2]);
    z           =   mxGetPr(prhs[3]);
    l           =   mxGetPr(prhs[4]);
    s           =   mxGetPr(prhs[5]);
    n_voxel     =   mxGetPr(prhs[6]);
    tbin1       =   mxGetScalar(prhs[7]);
    delta       =   mxGetScalar(prhs[8]);
    n_l         =   mxGetScalar(prhs[9]);
    n_s         =   mxGetScalar(prhs[10]);
    n_t         =   mxGetScalar(prhs[11]);
    z_wall      =   mxGetScalar(prhs[12]);

    plhs[0]     =   mxCreateDoubleMatrix(n_voxel[0]*n_voxel[1]*n_voxel[2], 1, mxREAL); 
    rho         =   mxGetPr(plhs[0]);

    for (int idx_l = 0; idx_l < n_l; idx_l++){
        if (idx_l%100 == 0){
            mexPrintf("%04d / %04d \n", idx_l, n_l);
            mexEvalString("drawnow;");
        }
        
        l_x = l[2*idx_l+0];
        l_y = l[2*idx_l+1];
        for (int idx_s = 0; idx_s < n_s; idx_s++){
            // mexPrintf("%d / %d \n", idx_s, n_s);
            // mexEvalString("drawnow;");
            s_x = s[2*idx_s+0];
            s_y = s[2*idx_s+1];

            idx_ls = n_t*(idx_s + n_s*idx_l);

            // find the first nonzero element for speedup
            firstNZ = -1;
            for (int idx_t = 0; idx_t < n_t; idx_t++){
                if(transient[idx_ls + idx_t] != 0){
                    firstNZ = idx_t;
                    break;
                }
            }

            if (firstNZ == -1){
                // mexPrintf("zero transient: %d / %d \n", idx_l, idx_s);
                // mexEvalString("drawnow;");
                continue; // if all zero, continue to next l,s pair
            }

            // find the last nonzero element for speedup
            lastNZ = n_t-1;
            for (int idx_t = n_t-1; idx_t >= 0; idx_t--){
                if(transient[idx_ls + idx_t] != 0){
                    lastNZ = idx_t;
                    lastz = z_wall - (idx_t*delta + tbin1)/2;
                    break;
                }
            }

            for (int idx_z = int(n_voxel[2])-1; idx_z >= 0; idx_z--){
                p_z = z[idx_z];
                if (p_z < lastz){
                    // mexPrintf("farther than last z: %f / %f \n", p_z, lastz);
                    // mexEvalString("drawnow;");
                    break; // if depth is far away than last depth, break to next l,s pair
                }
                if (p_z > z_wall){
                    continue;
                }

                for (int idx_y = 0; idx_y < int(n_voxel[1]); idx_y++){
                    p_y = y[idx_y];
                    for (int idx_x = 0; idx_x < int(n_voxel[0]); idx_x++){
                        p_x = x[idx_x];
                
                        r_l2 = pow(p_x - l_x, 2) + pow(p_y - l_y, 2) + pow(p_z - z_wall, 2);
                        r_l = sqrt(r_l2);
                        r_s2 = pow(p_x - s_x, 2) + pow(p_y - s_y, 2) + pow(p_z - z_wall, 2);
                        r_s = sqrt(r_s2);
                        r = r_l + r_s;
                        
                        idx_r = (r-tbin1)/delta;
                        idx_r_floor = floor(idx_r);
                        idx_r_ceil = idx_r_floor + 1;

                        idx_lst = idx_r_floor + idx_ls;
                        idx_zyx = idx_x + n_voxel[0]*idx_y + n_voxel[0]*n_voxel[1]*idx_z;
                        if (idx_r_floor >= firstNZ && idx_r_floor <= lastNZ){
                            rho[idx_zyx] += transient[idx_lst]*(double(idx_r_ceil) - idx_r); // (A>0)
                        }

                        if (idx_r_ceil >= firstNZ && idx_r_ceil <= lastNZ){
                            rho[idx_zyx] += transient[idx_lst+1]*(idx_r - double(idx_r_floor)); // (A>0)
                        }


                        // idx_r = int((r-tbin1)/delta + 0.5);
                        // if (idx_r < firstNZ || idx_r > lastNZ) {
                        //     continue; // if transient is 0, continue to next point
                        // }



                        // // if (idx_r < 0 || idx_r >= n_t) 
                        // //     continue; // if transient is 0, continue to next point
                        
                        // idx_lst = idx_r + idx_ls;

                        // idx_zyx = idx_x + n_voxel[0]*idx_y + n_voxel[0]*n_voxel[1]*idx_z;
                        // // if (transient[idx_lst] > 0.002){
                        // //     mexPrintf("transient!: (x, y, z) = (%f, %f, %f) / (l, s, t) = (%d, %d, %d) %f %f %f %f / %f\n", p_x, p_y, p_z, idx_l, idx_s, idx_r, transient[idx_lst]/0.002, r_l, r_s, r, lastz);
                        // //     mexEvalString("drawnow;");
                        // // }
                        // rho[idx_zyx] += transient[idx_lst];


                    }
                }
            }
        }
    }

	return;
}
