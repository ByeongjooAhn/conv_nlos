#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Check for proper number of arguments */

	if (nrhs != 13) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
				"mexComputeBackprojection requires 13 input arguments.");
	}

    // INPUT
    double  *x, *y, *l, *s, *n_voxel2;
    double  z_wall, x_target, y_target, z_target, z1, z_res;
    int     n_l, n_s;

    x_target    =   mxGetScalar(prhs[0]);
    y_target    =   mxGetScalar(prhs[1]);
    z_target    =   mxGetScalar(prhs[2]);
    x           =   mxGetPr(prhs[3]);
    y           =   mxGetPr(prhs[4]);
    l           =   mxGetPr(prhs[5]);
    s           =   mxGetPr(prhs[6]);
    z_wall      =   mxGetScalar(prhs[7]);
    n_voxel2    =   mxGetPr(prhs[8]);
    n_l         =   mxGetScalar(prhs[9]);
    n_s         =   mxGetScalar(prhs[10]);
    z1          =   mxGetScalar(prhs[11]);
    z_res       =   mxGetScalar(prhs[12]);


    // OUTPUT
    double  *kernel;

    plhs[0]     =   mxCreateDoubleMatrix(n_voxel2[0]*n_voxel2[1]*n_voxel2[2], 1, mxREAL); 
    kernel      =   mxGetPr(plhs[0]);

    // TEMP
    double p_x, p_y, l_x, l_y, s_x, s_y, L, S, A, D, z_ellipsoid, m, r_target_l, r_target_s, r_target_l2, r_target_s2, r_l2, r_s2, idx_z;
    int idx_zyx, idx_ls, idx_yx, idx_z_floor, idx_z_ceil;


    double *r_target    =   (double*) malloc(sizeof(double)*n_l*n_s);
    double *r_target2   =   (double*) malloc(sizeof(double)*n_l*n_s);
    double *m_target    =   (double*) malloc(sizeof(double)*n_l*n_s);

    // Compute target distance
    for (int idx_l = 0; idx_l < n_l; idx_l++){
        // mexPrintf("%d / %d \n", idx_l, n_l);
        // mexEvalString("drawnow;");
        l_x = l[2*idx_l+0];
        l_y = l[2*idx_l+1];
        for (int idx_s = 0; idx_s < n_s; idx_s++){
            s_x = s[2*idx_s+0];
            s_y = s[2*idx_s+1];
            idx_ls = idx_s + n_s*idx_l;

            r_target_l2 = pow(x_target - l_x, 2) + pow(y_target - l_y, 2) + pow(z_target - z_wall, 2);
            r_target_l = sqrt(r_target_l2);
            r_target_s2 = pow(x_target - s_x, 2) + pow(y_target - s_y, 2) + pow(z_target - z_wall, 2);
            r_target_s = sqrt(r_target_s2);

            r_target[idx_ls] = r_target_l + r_target_s;
            r_target2[idx_ls] = r_target[idx_ls]*r_target[idx_ls];
            if (r_target_l2 == 0 || r_target_s2 == 0){
                m_target[idx_ls] = 0;
            }
            else{
                m_target[idx_ls] = 1/(r_target_l2*r_target_s2);
            }
        }
    }


    // Compute kernel
    for (int idx_x = 0; idx_x < int(n_voxel2[0]); idx_x++){
        p_x = x[idx_x];
        for (int idx_y = 0; idx_y < int(n_voxel2[1]); idx_y++){
            p_y = y[idx_y];
        
            idx_yx = int(n_voxel2[2])*idx_y + int(n_voxel2[2])*int(n_voxel2[1])*idx_x;
            for (int idx_l = 0; idx_l < n_l; idx_l++){
                l_x = l[2*idx_l+0];
                l_y = l[2*idx_l+1];
                for (int idx_s = 0; idx_s < n_s; idx_s++){
                    s_x = s[2*idx_s+0];
                    s_y = s[2*idx_s+1];

                    idx_ls = idx_s + n_s*idx_l;
                   
                    L = pow(p_x - l_x, 2) + pow(p_y - l_y, 2);
                    S = pow(p_x - s_x, 2) + pow(p_y - s_y, 2);
                    A = (r_target[idx_ls]*r_target[idx_ls] - (L+S))/2;
                    D = (A*A - L*S)/(2*A + L + S);

                    if (D>0){
                    
                        z_ellipsoid = -sqrt(D) + z_wall;
                        idx_z = (z_ellipsoid-z1)/z_res;
                        idx_z_floor = floor(idx_z);
                        idx_z_ceil = idx_z_floor + 1;

                        idx_zyx = idx_z_floor + idx_yx;

                        m = m_target[idx_ls]; // A


                        // Anti-aliasing
                        if (idx_z_floor >= 0 && idx_z_floor < int(n_voxel2[2])){
                            kernel[idx_zyx] += m*(double(idx_z_ceil)-idx_z);// z->y->x order
                        }

                        if (idx_z_ceil >= 0 && idx_z_ceil < int(n_voxel2[2])){
                            kernel[idx_zyx+1] += m*(idx_z-double(idx_z_floor));// z->y->x order
                        }
                    }
                    // mexPrintf("%f / %f \n", double(idx_z_ceil)-idx_z, idx_z-double(idx_z_floor));
                    // mexEvalString("drawnow;");

                }
            }
        }
    }
    
    free(r_target);
    free(r_target2);
    free(m_target);

	return;
}
