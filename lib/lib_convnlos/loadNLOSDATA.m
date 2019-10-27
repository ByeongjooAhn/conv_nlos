function NLOSDATA = loadNLOSDATA(obj_name)
% loadNLOSDATA      Load light transient
%     Object is located at origin and LOS wall is at z = target_dist
% Input:
%     obj_name    Object name
% 
% Output:
%     DATA:
%     - obj_name                          :   Object name
%     - transient                         :   Light transient
%     - l                                 :   Illumination points
%     - s                                 :   Sensing points
%     - times                             :   Time bins
%     - is_confocal                       :   Scannning pattern
%     - target_dist                       :   Distance from LOS wall to center of target region center    
%     - target_size                       :   Width of target region
%     - delta                             :   ToF resolution (m)
%     - n_voxel                           :   Voxel resolution (for albedo)
%     - x / y / z                         :   Voxel position (for albedo)
%     - n_voxel_kernel                    :   Voxel resolution (for kernel)
%     - x_kernel / y_kernel / z_kernel    :   Voxel position (for kernel)
%     - n_sample_per_dim                  :   Kernel sampling resolution
%     - x_sample / y_sample / z_sample    :   Sampling location

%% Set parameters for each object
obj_list = {'S', 'USAF', 'soap', 'bunny', 'numbers', 'TX', '2019', 'toy'};

if any(strcmp(obj_name, obj_list))
    file_name = sprintf('./data/%s_transient.mat', obj_name);
    load(file_name);
else
    fprintf('Invalid obj_name\n');
end

end
