function kernel = FFTKernel(kernel, NLOSDATA, GPU_ENABLED)
% FFTKernel     Save FFT of kernel for efficient convolution

fprintf('FFTKernel\n');
tic;

% kernel basis
sum_kernel_basis = sum(kernel.basis(:)); % normalize kenrel for optimization
newsize = size(kernel_basis_1) + NLOSDATA.n_voxel - 1;

kernel_temp = kernel.basis/sum_kernel_basis;
if GPU_ENABLED
    kernel_temp = gpuArray(kernel_temp);
end
kernel.Fbasis = psf2otf(kernel_temp, newsize);

fprintf('Done! ');
toc;

end

