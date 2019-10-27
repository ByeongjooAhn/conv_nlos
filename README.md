# conv_nlos

This repsitory is an implementation of the method described in the following paper:

["Convolutional Approximations to the General Non-Line-of-Sight Imaging Operator"](http://imaging.cs.cmu.edu/conv_nlos/)\
[Byeongjoo Ahn](https://byeongjooahn.com), [Akshat Dave](https://ad74.blogs.rice.edu/), [Ashok Veeraraghavan](https://computationalimaging.rice.edu/team/ashok-veeraraghavan/), [Ioannis Gkioulekas](https://www.cs.cmu.edu/~igkioule/), and [Aswin C. Sankaranarayanan](https://users.ece.cmu.edu/~saswin/)\
IEEE/CVF International Conference on Computer Vision (ICCV), 2019    (Oral Presentation)



## Run Code

1. Run `lib/lib_convnlos/mex/complie_mex.m` to complie mex codes.
2. Run `demo.m` script. It includes Light Cone Transform (LCT), filtered backprojection, our method without priors, and our method with priors for comparison. 
3. Results are stored in `results/[obj_name]/`. It includes runtime for each algorithm, and reconstructed NLOS albedo and its 2D maximum intensity projection in three different colormaps.



## Inputs

Inputs are stored in `data/[obj_name]_transient.mat` .  It includes  `NLOSDATA` structure consisting of the following parameters.

| Parameter                      | Description                                                  |
| ------------------------------ | ------------------------------------------------------------ |
| obj_name                       | Object name.                                                 |
| transient                      | Light transient [n_l * n_s * n_t]. Transients need to be rectified such that the direct component starts at t = 0. |
| l                              | Illumination points [n_l * 3].                               |
| s                              | Sensing points [n_s *3].                                     |
| times                          | Time bins [1* n_t].                                          |
| is_confocal                    | Scanning pattern (1 for confocal, 0 for 5D).                 |
| target_dist                    | Distance from LOS wall to object center (m). Object is located at origin, and LOS wall is at z = target_dist. |
| delta                          | ToF resolution (m).                                          |
| n_voxel                        | Voxel resolution for NLOS albedo [n_x * n_y * n_z].          |
| x / y / z                      | Voxel position for NLOS albedo [1 * n_x] / [1 * n_y] / [1 * n_z]. |
| n_voxel_kernel                 | Voxel resolution for kernel.                                 |
| x_kernel / y_kernel / z_kernel | Voxel position for kernel.                                   |
| n_sample_per_dim               | Kernel sampling resolution. We use n_sample_per_dim = 3, which means the kernels are sampled from 3 * 3 * 3 = 27 points to make a kernel basis. |
| x_sample / y_sample / z_sample | Kernel sampling location.                                    |

