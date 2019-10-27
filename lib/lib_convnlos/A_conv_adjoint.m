function rho_out = A_conv_adjoint(rho_in, Fkernel_basis, n_voxel)
    
    rho_in = reshape(rho_in, n_voxel);
    frho_in = fftn(rho_in, size(Fkernel_basis(:,:,:,1)));

    krho_temp = ifftn(frho_in.*conj(Fkernel_basis));
    rho_out = krho_temp(1:n_voxel(1), 1:n_voxel(2), 1:n_voxel(3));
    
    rho_out = real(rho_out(:));
end
   