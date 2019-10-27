function NLOSDATA = loadNLOSDATA_REVEAL(obj_name)
% loadNLOSDATA      Load light transient
% Input:
%     obj_name    Object name
% 
% Output:
%     DATA:
%     - obj_name                          :   Object name
%     - transient
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
%     - SNR                               :   SNR for LCT
%     - bias                              :   ToF bias for correcting delay

%     - x_coeff, y_coeff, z_coeff         :   Size of target region / Size of wall
%     - K                                 :   Downsampling
%     - n_point_per_dim                   :   Voxel resolution per axis
%     - n_sample_per_dim                  :   Kernel sampling resolution
%     - x_sample / y_sample / z_sample    :   Sampling location

%% Set parameters for each object
% obj_list = {'4', 'bowl_front', '7spad', '4_10s', '7spad_2'};
% 
% if any(strcmp(obj_name, obj_list))
%     file_name = sprintf('./data_REVEAL/%s_transient.mat', obj_name);
%     load(file_name);
% else
%     fprintf('Invalid obj_name\n');
% end
% 
% 
% obj_name = '4_10s';
% width = 1.5;
% tof_resolution = 8;
% nPoints = 128;
% K = 3; % 8 -> 64 ps
% 
% NLOSDATA = loadFromREVEAL(obj_name, mat, tof, nPoints, width, tof_resolution, K);
% save(sprintf('./data_REVEAL/%s_transient.mat', obj_name), 'NLOSDATA');
% 
% % 
num_spad = 1;
dataset = cell(num_spad, 1);
% for i = 1:num_spad
%     dataset{i} = load(sprintf('data_REVEAL_raw/7spad/data_spad_%d.mat', i));
% end

% % %
spad_id = 1;
dataset{1} = load(sprintf('data_REVEAL_raw/7spad/data_spad_%d.mat', spad_id));
% % % 

obj_name = sprintf('7spad_2_%d', spad_id);
K = 4; % K=4: 4 -> 64ps
NLOSDATA = loadFrom7SPAD(obj_name, dataset, K, num_spad);
save(sprintf('./data_REVEAL/%s_transient.mat', obj_name), 'NLOSDATA');

% [rho, time, NLOSDATA] = computeBPFrom7SPAD(obj_name, dataset, K);

end

% code snippet for saving



