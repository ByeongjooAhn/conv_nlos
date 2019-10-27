function saveAlbedo_cmap(rho, METHOD, DIRNAME)

load('lib/util/cmap.mat');

rho = real(rho);

xy = permute(squeeze(max(rho(:,end:-1:1,:), [], 3)), [2 1]);
imwrite(ind2rgb(im2uint8(mat2gray(xy)), colormap(cmap)), sprintf('%s/xy_%s.png', DIRNAME, METHOD));

xz = permute(squeeze(max(rho, [], 2)), [2 1]);
imwrite(ind2rgb(im2uint8(mat2gray(xz)), colormap(cmap)), sprintf('%s/xz_%s.png', DIRNAME, METHOD));

yz = permute(squeeze(max(rho(:,end:-1:1,:), [], 1)), [1 2]);
imwrite(ind2rgb(im2uint8(mat2gray(yz)), colormap(cmap)), sprintf('%s/yz_%s.png', DIRNAME, METHOD));

end
