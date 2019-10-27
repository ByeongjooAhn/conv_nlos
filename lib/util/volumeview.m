function y = volumeview(x)
    
    width = 1;
    [val, idx] = max(x, [], 3);

    mip = zeros(size(x));
    for idx_x = 1:size(idx,1)
        for idx_y = 1:size(idx,2)
            idx_temp = idx(idx_x, idx_y);
            idx_temp_m = max(idx_temp - width, 1);
            idx_temp_p = min(idx_temp + width, size(x, 3));
            mip(idx_x, idx_y, idx_temp_m:idx_temp_p) = val(idx_x, idx_y);
        end
    end
    
    y = mip/max(mip(:));
    
end