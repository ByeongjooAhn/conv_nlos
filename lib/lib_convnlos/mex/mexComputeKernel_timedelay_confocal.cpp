#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Check for proper number of arguments */

	if (nrhs != 11) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
				"mexComputeBackprojection requires 11 input arguments.");
	}

    // INPUT
    double  *x, *y, *l, *n_voxel2;
    double  z_wall, x_target, y_target, z_target, z1, z_res;
    int     POINTS;

    x_target    =   mxGetScalar(prhs[0]);
    y_target    =   mxGetScalar(prhs[1]);
    z_target    =   mxGetScalar(prhs[2]);
    x           =   mxGetPr(prhs[3]);
    y           =   mxGetPr(prhs[4]);
    l           =   mxGetPr(prhs[5]);
    z_wall      =   mxGetScalar(prhs[6]);
    n_voxel2    =   mxGetPr(prhs[7]);
    POINTS      =   mxGetScalar(prhs[8]);
    z1          =   mxGetScalar(prhs[9]);
    z_res       =   mxGetScalar(prhs[10]);


    // OUTPUT
    double  *kernel;

    plhs[0]     =   mxCreateDoubleMatrix(n_voxel2[0]*n_voxel2[1]*n_voxel2[2], 1, mxREAL);
    kernel      =   mxGetPr(plhs[0]);

    // TEMP
    double p_x, p_y, l_x, l_y, s_x, s_y, L, A, D, z_ellipsoid, m, r_target_l, r_target_l2, r_l2, idx_z;
    int idx_zyx, idx_yx, idx_z_floor, idx_z_ceil;


    double *r_target    =   (double*) malloc(sizeof(double)*POINTS);
    double *r_target2   =   (double*) malloc(sizeof(double)*POINTS);
    double *m_target    =   (double*) malloc(sizeof(double)*POINTS);

    // Compute target distance
    for (int idx_l = 0; idx_l < POINTS; idx_l++){
        // mexPrintf("%d / %d \n", idx_l, POINTS);
        // mexEvalString("drawnow;");
        l_x = l[2*idx_l+0];
        l_y = l[2*idx_l+1];

        r_target_l2 = pow(x_target - l_x, 2) + pow(y_target - l_y, 2) + pow(z_target - z_wall, 2);
        r_target_l = sqrt(r_target_l2);

        r_target[idx_l] = 2*r_target_l;
        r_target2[idx_l] = r_target[idx_l]*r_target[idx_l];
        m_target[idx_l] = 1/(r_target_l2*r_target_l2);

    }


    // Compute kernel
    for (int idx_y = 0; idx_y < int(n_voxel2[1]); idx_y++){
        p_y = y[idx_y];
        for (int idx_x = 0; idx_x < int(n_voxel2[0]); idx_x++){
            p_x = x[idx_x];

            idx_yx = int(n_voxel2[2])*idx_x + int(n_voxel2[2])*int(n_voxel2[1])*idx_y;
            for (int idx_l = 0; idx_l < POINTS; idx_l++){

                l_x = l[2*idx_l+0];
                l_y = l[2*idx_l+1];

                L = pow(p_x - l_x, 2) + pow(p_y - l_y, 2);
                D = (pow(r_target2[idx_l] - L - L, 2) - 4*L*L)/(4*r_target2[idx_l]);

                if (D>0){
                    
                    z_ellipsoid = -sqrt(D) + z_wall;
                    idx_z = (z_ellipsoid-z1)/z_res;
                    idx_z_floor = floor(idx_z);
                    idx_z_ceil = idx_z_floor + 1;

                    idx_zyx = idx_z_floor + idx_yx;

                    // m = 1/(r_l2*r_s2); // A
                    m = m_target[idx_l]; // A>0

                    // Anti-aliasing
                    if (idx_z_floor >= 0 && idx_z_floor < int(n_voxel2[2])){
                        kernel[idx_zyx] += m*(double(idx_z_ceil)-idx_z);// z->y->x order
                    }

                    if (idx_z_ceil >= 0 && idx_z_ceil < int(n_voxel2[2])){
                        kernel[idx_zyx+1] += m*(idx_z-double(idx_z_floor));// z->y->x order
                    }
                }
            }
        }
    }

    free(r_target);
    free(r_target2);
    free(m_target);

	return;

    
}
