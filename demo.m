%% "Convolutional Approximations to the General NLOS Imaging operator"
% Compatible with MATLAB 2016a or higher (lower versions have not been tested)
% Run ./lib/lib_convnlos/mex/compile_mex.m before running this code
% Email bahn@cmu.edu if you have any issues using this code

clear
close all
restoredefaultpath;
addpath(genpath('./lib'));

%% SETUP
% Use GPU if possible (only for convolution)
GPU_ENABLED = logical(gpuDeviceCount);
fprintf('GPU_ENABLED: %d\n', GPU_ENABLED);

%% LOAD DATA
% Select an object from obj_list
obj_list = {'S', 'USAF', 'soap', 'bunny', 'numbers', 'TX', '2019', 'toy'};
obj_name = obj_list{7};
NLOSDATA = loadNLOSDATA(obj_name);

%% RUN BASELINES
% LCT
[rho.LCT, time.LCT] = computeLCT(NLOSDATA, GPU_ENABLED);

% Filtered backprojection
is_load_BP = false;
[rho.BP, time.BP] = computeBP(NLOSDATA, is_load_BP); % alternative BP
[rho.FBP, time.FBP] = computeFBP(rho.BP);

%% RUN OURS
% Make kernel
is_load_kernel = true;
kernel = makeKernel(NLOSDATA, GPU_ENABLED, is_load_kernel);

% Ours w/o priors
[rho.Gram_noprior, time.Gram_noprior] = computeOurs_noprior(rho.BP, kernel, GPU_ENABLED);

% Ours (w/ priors) - change weight_l1 15e-2 ~ 20e-2
params_admm.mu = 1;
params_admm.iter = 40;
params_admm.nonneg = true;
params_admm.weight_l1 = 20e-2;
params_admm.weight_TV = 1e-2;
rho_init = rho.Gram_noprior;
[rho.Gram, time.Gram] = computeOurs(rho.BP, rho_init, kernel, params_admm, NLOSDATA, GPU_ENABLED);

figure; imslice(rho.Gram);

%% SAVE RESULTS
saveResults(rho, time, NLOSDATA, kernel);
