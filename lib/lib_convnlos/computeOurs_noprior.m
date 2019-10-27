function [rho, time] = computeOurs_noprior(rho_bp, kernel, GPU_ENABLED)
% computeOurs_noprior   Compute ours without prior

fprintf('[computeOurs_noprior]\n');
tic;

% Set function
n_voxel = size(rho_bp);
Afun = @(rho_in) A_conv(rho_in(:), kernel.Fbasis, n_voxel);

if GPU_ENABLED
     rho_bp = gpuArray(rho_bp);
end

% Run optimization
params.tol = 1e-2;
params.iter = 100;
rho = cgs(Afun, rho_bp(:), params.tol, params.iter);

rho = reshape(rho, size(rho_bp));
rho = real(rho);
rho(rho<0) = 0;

% Edge
rho(:,:,1) = 0;
rho(:,:,end) = 0;
rho = rho./max(rho(:));
rho = gather(rho);

time = toc;
fprintf('Done! ');
toc;

end

