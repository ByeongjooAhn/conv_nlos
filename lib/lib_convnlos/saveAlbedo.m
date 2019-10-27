function saveAlbedo(rho, METHOD, DIRNAME)

rho = real(rho);

xy = permute(squeeze(max(rho(:,end:-1:1,:), [], 3)), [2 1]);
imwrite(ind2rgb(im2uint8(mat2gray(xy)), parula(256)), sprintf('%s/xy_%s.png', DIRNAME, METHOD));

xz = permute(squeeze(max(rho, [], 2)), [2 1]);
imwrite(ind2rgb(im2uint8(mat2gray(xz)), parula(256)), sprintf('%s/xz_%s.png', DIRNAME, METHOD));

yz = permute(squeeze(max(rho(:,end:-1:1,:), [], 1)), [1 2]);
imwrite(ind2rgb(im2uint8(mat2gray(yz)), parula(256)), sprintf('%s/yz_%s.png', DIRNAME, METHOD));

end
