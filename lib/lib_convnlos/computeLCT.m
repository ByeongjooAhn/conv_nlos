function [rho, time] = computeLCT(NLOSDATA, GPU_ENABLED)
% computeLCT    Compute light cone transform
%     Modified input format from the original LCT code

fprintf('[computeLCT] ');

SNR = 8e-1;

if NLOSDATA.is_confocal
    transient_ = NLOSDATA.transient;
else % extract confocal transient if input is 5D
    transient_ = NLOSDATA.transient_confocal;
end

delta = NLOSDATA.delta;
times = delta:delta:NLOSDATA.times(end);
M = length(times);

% Pad zero
[N, ~, M_] = size(transient_);
if M <= M_
    transient_confocal = transient_;
    M = M_;
    times = linspace(delta, NLOSDATA.times(end), M);
else
    transient_confocal = zeros(N, N, M);
    transient_confocal(:,:,end-M_+1:end) = transient_;
end
data = permute(transient_confocal, [3 2 1]);
[mtx,mtxi] = resamplingOperator(M);

% Set parameters
range = times(end);
width = (max(NLOSDATA.l(:,1)) - min(NLOSDATA.l(:,1)))/2;

% Define kernel
psf = definePSF(N, M, width/range);

tic;
if GPU_ENABLED
    data = gpuArray(data);
    psf = gpuArray(psf);
end
fpsf = fftn(psf);
invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./SNR);

% Define transform operators
grid_z = repmat(times',[1 N N]);
data = data.*(grid_z).^4; % (isDiffuse) (c1/sqrt(x) included in mtx)

% Step 2: Resample time axis and pad result
tdata = zeros(2*M-1,2*N-1,2*N-1);
if GPU_ENABLED
    tdata = gpuArray(tdata);
end
tdata(1:M,1:N,1:N) = reshape(mtx*data(:,:),[M N N]);

% Step 3: Convolve with inverse filter and unpad result
tvol = ifftn(fftn(tdata).*invpsf);
tvol = tvol(1:M,1:N,1:N);

% Step 4: Resample depth axis and clamp results
vol  = reshape(mtxi*tvol(:,:),[M N N]);
vol  = max(real(vol),0);

rho = permute(vol, [3 2 1]);

% Crop for same input size
z_min = min(NLOSDATA.z);
z_max = max(NLOSDATA.z);
start_idx = round(((NLOSDATA.target_dist - z_max)*2)  / NLOSDATA.delta);
start_idx = max(start_idx, 1);
end_idx = round(((NLOSDATA.target_dist - z_min)*2) / NLOSDATA.delta);
rho(:,:,end+1:end_idx) = 0;

rho = rho(:,:,start_idx:end_idx);
rho = rho(:,:,end:-1:1); % align with our format (far to near)
rho = rho/max(rho(:));
rho = gather(rho);

time = toc; 
fprintf('Done! ');
toc;

end


function psf = definePSF(N, M, slope)
% Local function to compute NLOS blur kernel

x = linspace(-1,1,2*N-1);
y = linspace(-1,1,2*N-1);
z = linspace(0,2,2*M-1);
[grid_z,grid_y,grid_x] = ndgrid(z,y,x);

% Define PSF
psf = abs(((4*slope).^2).*(grid_x.^2 + grid_y.^2) - grid_z);
psf = double(psf == repmat(min(psf,[],1),[2*M-1 1 1]));
psf = psf./norm(psf(:));
psf = circshift(psf,[0 N N]);

end

% function [mtx,mtxi] = resamplingOperator(M) % (original)
% % Local function that defines resampling operators
% 
% mtx = sparse([],[],[],M.^2,M,M.^2);
% 
% x = 1:M.^2;
% mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
% mtx  = spdiags(1./sqrt(x)',0,M.^2,M.^2)*mtx;
% mtxi = mtx';
% 
% K = log(M)./log(2);
% for k = 1:round(K)
%     mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
%     mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
% end
% end

function [mtx,mtxi] = resamplingOperator(M) % (updated by Byeongjoo)
% Local function that defines resampling operators

mtx = sparse([],[],[],M^2,M,M^2);

x = 1:M^2;
mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
mtx  = spdiags(1./sqrt(x)',0,M^2,M^2)*mtx;

K = kron(speye(M), ones(1, M)); 

mtx = K*mtx;
mtxi = mtx';

end
