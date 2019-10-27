function [rho, time] = computeBP(NLOSDATA, is_load)
% ComputeBP     Compute backprojection (alternative)
%     Alternative BP corresponding to (A>0)' operator

fprintf('[computeBP]\n');

bp_dir = sprintf('data/temp/%s/bp.mat', NLOSDATA.obj_name);
mkdir(sprintf('data/temp/%s', NLOSDATA.obj_name))
if is_load
    if exist(bp_dir, 'file') == 2
        load(bp_dir);

        fprintf('Loaded!\n');

        return;
    end
end
tic;

% Set target space
x = NLOSDATA.x;
y = NLOSDATA.y;
z = NLOSDATA.z;

l = NLOSDATA.l;
n_t = length(NLOSDATA.times);
tbin1 = NLOSDATA.times(1);
z_wall = l(1,3);

transient = permute(NLOSDATA.transient, [3 2 1]);
transient = transient(:);

% Backprojection (Input for mex must be double!)
if NLOSDATA.is_confocal
    rho = mexComputeBackprojection_timedelay_confocal(transient, x, y, z, l(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, size(NLOSDATA.l,1), n_t, z_wall);
    rho = reshape(rho, NLOSDATA.n_voxel);
else
    s = NLOSDATA.s;
    rho = mexComputeBackprojection_timedelay(transient, x, y, z, l(:,1:2)', s(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, size(NLOSDATA.l,1), size(NLOSDATA.s,1), n_t, z_wall);
    rho = reshape(rho, NLOSDATA.n_voxel);
end

rho = rho/max(rho(:));

time = toc;
fprintf('Done! ');
toc;

save(bp_dir, 'rho', 'time');

end
