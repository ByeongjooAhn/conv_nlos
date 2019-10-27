function rho = computeBP_iter(transient_in, NLOSDATA)
% computeBP_iter Compute backprojection (alternative) for iteration
%     Slightly different output format for iteration
    
% Set target space
x = NLOSDATA.x;
y = NLOSDATA.y;
z = NLOSDATA.z;

l = NLOSDATA.l;
n_t = length(NLOSDATA.times);
tbin1 = NLOSDATA.times(1);
z_wall = l(1,3);

transient = permute(transient_in, [3 2 1]);
transient = transient(:);

% Input should be double!
if NLOSDATA.is_confocal
    rho = mexComputeBackprojection_timedelay_confocal(transient, x, y, z, l(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, NLOSDATA.n_point_per_dim, n_t, z_wall);
else
    s = NLOSDATA.s;
    rho = mexComputeBackprojection_timedelay(transient, x, y, z, l(:,1:2)', s(:,1:2)', NLOSDATA.n_voxel, tbin1, NLOSDATA.delta, size(NLOSDATA.l,1), size(NLOSDATA.s,1), n_t, z_wall);
end

rho = reshape(rho, NLOSDATA.n_voxel);

end
