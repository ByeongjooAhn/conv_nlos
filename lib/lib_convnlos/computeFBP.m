function [rho, time] = computeFBP(rho_bp)
% computeFBP    Compute filtered backprojection (Laplacian)

fprintf('[computeFBP] ');
tic;

dz = zeros(3,3,3);
dz(2,2,1) = -1;
dz(2,2,2) = 2;
dz(2,2,3) = -1;

newsize = size(rho_bp) + size(dz) - 1;
d3FT = psf2otf(dz, newsize);

K3 = @(x) real(ifftn(d3FT .* fftn(x, newsize)));

rho = K3(rho_bp);
rho = rho(1:size(rho_bp,1), 1:size(rho_bp,2), 1:size(rho_bp,3)); 

% edge case
rho(:,:,1) = 0;
rho(:,:,end) = 0;

rho(rho<0) = 0;
rho = rho/max(rho(:));

time = toc;
fprintf('Done! ');
toc;

end