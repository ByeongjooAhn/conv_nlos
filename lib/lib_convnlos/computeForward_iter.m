function [transient] = computeForward_iter(rho, NLOSDATA)
% computeForward_iter Compute forward operation for iteration

% Set target space
x = NLOSDATA.x;
y = NLOSDATA.y;
z = NLOSDATA.z;

l = NLOSDATA.l;
n_t = length(NLOSDATA.times);
tbin1 = NLOSDATA.times(1);
z_wall = l(1,3);

rho = rho(:);

% Input should be double!
if NLOSDATA.is_confocal
    transient = mexComputeForward_timedelay_confocal(rho, x, y, z, l(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, NLOSDATA.n_point_per_dim, n_t, z_wall);
    transient = reshape(transient, [n_t NLOSDATA.n_point_per_dim, NLOSDATA.n_point_per_dim]);
    transient = permute(transient, [3 2 1]);
else
    s = NLOSDATA.s;
    transient = mexComputeForward_timedelay(rho, x, y, z, l(:,1:2)', s(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, size(NLOSDATA.l,1), size(NLOSDATA.s,1), n_t, z_wall);
    transient = reshape(transient, [n_t size(l,1), size(s,1)]);
    transient = permute(transient, [3 2 1]);
end

end
