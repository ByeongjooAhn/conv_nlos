function [rho, time] = computeFull(rho_bp, NLOSDATA)
% computeFull   Compute full linear reconstruction 

fprintf('computeFull ');
tic;

params.tol = 1e-2;
params.iter = 40;

Afun = @(rho_in) computeAfunFull(rho_in(:), NLOSDATA);

rho = cgs(Afun, rho_bp(:), params.tol, params.iter);
rho = reshape(rho, size(rho_bp));
rho = real(rho);

% Edge
rho(:,:,1) = 0;
rho(:,:,end) = 0;

rho = rho./max(rho(:));

time = toc;
fprintf('Done! ');
toc;

end


function rho_out = computeAfunFull(rho_in, NLOSDATA)
    transient = computeForward_iter(rho_in, NLOSDATA);
    rho_out = computeBP_iter(transient, NLOSDATA);
    rho_out = rho_out(:);
    fprintf('computeAfunFull Done!\n');
end

