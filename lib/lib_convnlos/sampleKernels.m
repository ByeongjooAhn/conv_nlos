function kernel = sampleKernels(NLOSDATA)
% sampleKernels     Sample kernels to make a rank-1 kernel 

fprintf('sampleKernels\n');
tic;
n_sample_per_dim = NLOSDATA.n_sample_per_dim;

x_sample = NLOSDATA.x_sample;
y_sample = NLOSDATA.y_sample;
z_sample = NLOSDATA.z_sample;

z_res = NLOSDATA.z(2) - NLOSDATA.z(1);
z_wall = NLOSDATA.l(1,3);

%% SAMPLE KERNELS
kernel_sampled = zeros([prod(NLOSDATA.n_voxel_kernel) n_sample_per_dim^3]);
for idx_x = 1:n_sample_per_dim
    fprintf('Sampling kernels... %d/%d\n', idx_x, n_sample_per_dim);
    for idx_y = 1:n_sample_per_dim

        x_target = x_sample(idx_x);
        y_target = y_sample(idx_y);
        x = x_target + NLOSDATA.x_kernel;
        y = y_target + NLOSDATA.y_kernel;
        
        for idx_z = 1:n_sample_per_dim
            z_target = z_sample(idx_z);
            z = z_target + NLOSDATA.z_kernel;
                        
            if NLOSDATA.is_confocal
                kernel_temp = mexComputeKernel_timedelay_confocal(x_target, y_target, z_target, x, y, NLOSDATA.l(:,1:2)', z_wall, NLOSDATA.n_voxel_kernel, NLOSDATA.n_point_per_dim, z(1), z_res);
                kernel_temp = kernel_temp*(NLOSDATA.target_dist-z_target)^2;
            else
%                 kernel_temp = mexComputeKernel_timedelay(x_target, y_target, z_target, x, y, NLOSDATA.l(:,1:2)', NLOSDATA.s(:,1:2)', z_wall, NLOSDATA.n_voxel_kernel, size(NLOSDATA.l,1), size(NLOSDATA.s,1), z(1), z_res);
                kernel_temp = mexComputeKernel_timedelay_aa(x_target, y_target, z_target, x, y, NLOSDATA.l(:,1:2)', NLOSDATA.s(:,1:2)', z_wall, NLOSDATA.n_voxel_kernel, size(NLOSDATA.l,1), size(NLOSDATA.s,1), z(1), z_res);
            end
            
            kernel_temp = reshape(kernel_temp, NLOSDATA.n_voxel_kernel(end:-1:1));
            kernel_temp = permute(kernel_temp, [3 2 1]);
            kernel_temp = kernel_temp(:);
            
            idx_zxy = idx_y + n_sample_per_dim*(idx_x-1) + n_sample_per_dim^2*(idx_z-1);
            kernel_sampled(:, idx_zxy) = kernel_temp;
        end
    end
end

%% MAKE BASIS
fprintf('Make basis...\n');

% Normalize kernels_sampled
kernels_sampled_norm = vecnorm(kernel_sampled, Inf); % Inf norm ~ z^-4
kernels_sampled_normalized = kernel_sampled./kernels_sampled_norm;
kernels_sampled_normalized(isnan(kernels_sampled_normalized)) = 0;

% SVD
[U, ~, ~] = svd(kernels_sampled_normalized, 0);
kernel_sampled = reshape(kernels_sampled_normalized, [NLOSDATA.n_voxel_kernel n_sample_per_dim^3]);

% Make kernel positive
sign_U = sign(sum(U(:,1)));
U = sign_U*U;

kernel_basis = reshape(U, [NLOSDATA.n_voxel_kernel n_sample_per_dim^3]);

%% EXTRACT KERNELS
kernel.basis = kernel_basis(:,:,:,1);
kernel.sampled = kernel_sampled;
kernel.singular_values = singular_values;

fprintf('Done! ');
toc;

end
