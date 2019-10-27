function imslice(x)

subplot 224; imagesc(permute(squeeze(max(x(:,end:-1:1,:), [], 3)), [2 1])); 
axis image; xlabel('x'); ylabel('y'); % colorbar;
title('Front view XY (MIP)');

subplot 222; imagesc(permute(squeeze(max(x, [], 2)), [2 1])); 
axis image; xlabel('x'); ylabel('z'); % colorbar;
title('Top view XZ (MIP)');

subplot 223; imagesc(permute(squeeze(max(x(:,end:-1:1,:), [], 1)), [1 2])); 
axis image; xlabel('z'); ylabel('y'); % colorbar;
title('Left view ZY (MIP)');

x = volumeview(x);
x_permute = permute(x, [1 3 2]);
subplot 221; vol3d('CData', x_permute); 
axis image; xlabel('z'); ylabel('x'); zlabel('y'); colorbar;

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])
set(gca, 'Xcolor', 'none');
set(gca, 'Ycolor', 'none');
set(gca, 'Zcolor', 'none');

view(135,35);

set(gca, 'color', 'k');
colormap gray;

end