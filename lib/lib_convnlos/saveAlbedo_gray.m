function saveAlbedo_gray(rho, METHOD, DIRNAME)

rho = real(rho);

xy = permute(squeeze(max(rho(:,end:-1:1,:), [], 3)), [2 1]);
xy = xy/max(xy(:));
imwrite(xy, sprintf('%s/xy_%s.png', DIRNAME, METHOD));

xz = permute(squeeze(max(rho, [], 2)), [2 1]);
xz = xz/max(xz(:));
imwrite(xz, sprintf('%s/xz_%s.png', DIRNAME, METHOD));

yz = permute(squeeze(max(rho(:,end:-1:1,:), [], 1)), [1 2]);
yz = yz/max(yz(:));
imwrite(yz, sprintf('%s/yz_%s.png', DIRNAME, METHOD));

end
