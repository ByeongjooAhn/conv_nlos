function [rho, time] = computeOurs(rho_bp, rho_init, kernel, params, NLOSDATA, GPU_ENABLED)
% computeOurs   Compute ours with priors
%     Incorporates positivity, sparsity, and total variation
%     Modified from linearized admm of LCT code

fprintf('[computeOurs]\n');
tic;

n_voxel = size(rho_bp);

% quadratic operator and transpose
shrinkage = @(a,kappa) repmat(max(0, 1-kappa./sqrt(a(:,:,:,1).^2 + a(:,:,:,2).^2 + a(:,:,:,3).^2)),1,1,1,3).*a;

max_iters = params.iter;

% Initialize ADMM variables
u1 = zeros(n_voxel);
u2 = zeros(n_voxel);
u3 = zeros(n_voxel);
u4 = zeros(n_voxel);
u5 = zeros(n_voxel);
u6 = zeros(n_voxel);
 
% Define the A operator 
K1 = @(x) reshape(A_conv(x(:), kernel.Fbasis, n_voxel), n_voxel);

% Define the A_adjoint operator 
K1T = @(x) reshape(A_conv_adjoint(x(:), kernel.Fbasis, n_voxel), n_voxel);

% gradient kernels
d2 = [0 0 0; 0 1 -1; 0 0 0];
d2 = padarray(d2, [0,0,1]);
d1 = [0 0 0; 0 1 0; 0 -1 0];
d1 = padarray(d1, [0,0,1]);
d3 = zeros(3,3,3);
d3(2,2,2) = 1; d3(2,2,3) = -1;

% operator functions
p2o = @(x) psf2otf(x, n_voxel);
d1FT = p2o(d1);
d2FT = p2o(d2);
d3FT = p2o(d3);

K2 = @(x) x;
K3 = @(x) real(ifftn(d1FT .* fftn(x))); %debug
K4 = @(x) real(ifftn(d2FT .* fftn(x)));
K5 = @(x) real(ifftn(d3FT .* fftn(x)));
K6 = @(x) x;

K2T = @(x) x;
K3T = @(x) real(ifftn(conj(d1FT) .* fftn(x)));
K4T = @(x) real(ifftn(conj(d2FT) .* fftn(x)));
K5T = @(x) real(ifftn(conj(d3FT) .* fftn(x)));
K6T = @(x) x;

% the full operator and transpose
K = @(x) [K1(x); ...
    K2(x); ...
    K3(x); ...
    K4(x);...
    K5(x);...
    K6(x)];

KT = @(x) K1T(x(1:end/6,:,:))...
    + K2T(x(end/6+1:2*end/6,:,:)) ...
    + K3T(x(2*end/6+1:3*end/6,:,:)) ...
    + K4T(x(3*end/6+1:4*end/6,:,:))...
    + K5T(x(4*end/6+1:5*end/6,:,:))...
    + K6T(x(5*end/6+1:6*end/6,:,:));

if GPU_ENABLED
    K1_CPU = @(x) reshape(A_conv(x(:), gather(kernel.Fbasis), n_voxel), n_voxel);
    K1T_CPU = @(x) reshape(A_conv_adjoint(x(:), gather(kernel.Fbasis), n_voxel), n_voxel);

    K_CPU = @(x) [K1_CPU(x); ...
        K2(x); ...
        K3(x); ...
        K4(x);...
        K5(x);...
        K6(x)];
    KT_CPU = @(x) K1T_CPU(x(1:end/6,:,:))...
        + K2T(x(end/6+1:2*end/6,:,:)) ...
        + K3T(x(2*end/6+1:3*end/6,:,:)) ...
        + K4T(x(3*end/6+1:4*end/6,:,:))...
        + K5T(x(4*end/6+1:5*end/6,:,:))...
        + K6T(x(5*end/6+1:6*end/6,:,:));
else
    K_CPU = K;
    KT_CPU = KT;
end

% Compute and save operator norm (slow)
norm_dir = sprintf('data/%s_operator_norm.mat', NLOSDATA.obj_name);
if exist(norm_dir, 'file') == 2
    load(norm_dir);
else
%     fprintf('computing operator norm (only for first run) ... ');
%     [Knorm, Ksmall] = compute_operator_norm(K_CPU, KT_CPU, n_voxel);
%     save(norm_dir, 'Knorm', 'Ksmall');
%     fprintf('Done!\n');
    
    Knorm = 3.74; % empirical
end

% fprintf('Operator eigenvalue ratio: %.02f (%.02f/%.02f)\n', Knorm/Ksmall, Knorm, Ksmall);

b = rho_bp;
x = rho_init;

if GPU_ENABLED
    x = gpuArray(x);
end

% Run linearized ADMM
mu = params.mu;
nu = mu * Knorm^2;
lambda_l1 = mu * params.weight_l1;
lambda_TV = mu * params.weight_TV;

for ii = 1:max_iters
   fprintf('%02d / %02d\n', ii, max_iters);

   % z1 update -- quadratic term
   v = K1(x) + u1; 
   z1 = (b + mu*v)/(1+mu);
   
   % z2 update -- non-negativity
   v = x + u2;
   if params.nonneg
        z2 = max(0,v); 
   else
        z2 = v;
   end

   % z3/z4/z5 update -- TV regularizer
   v = cat(4, K3(x), K4(x), K5(x)) + cat(4, u3, u4, u5);
   kappa = lambda_TV/mu;
   shrunk = shrinkage(v, kappa);
   z3 = squeeze(shrunk(:,:,:,1));
   z4 = squeeze(shrunk(:,:,:,2));
   z5 = squeeze(shrunk(:,:,:,3));

   % z6 update -- Sparsity penalty
   v = x + u6;
   kappa = lambda_l1/mu;
   z6 = max(v - kappa,0) - max(-v - kappa,0); 

   % u update
   u1 = u1 + K1(x) - z1;
   u2 = u2 + x - z2;
   u3 = u3 + K3(x) - z3;
   u4 = u4 + K4(x) - z4;
   u5 = u5 + K5(x) - z5;
   u6 = u6 + K6(x) - z6;

   % x update
   v1 = z1 - u1;
   v2 = z2 - u2;
   v3 = z3 - u3;
   v4 = z4 - u4;
   v5 = z5 - u5;
   v6 = z6 - u6;
   x = x - (mu/nu) * (KT(K(x)) - KT([v1; v2; v3; v4; v5; v6]));
end

rho = real(x);
rho = rho./max(rho(:));
rho = gather(rho);

time = toc;
fprintf('Done! ');
toc;

end