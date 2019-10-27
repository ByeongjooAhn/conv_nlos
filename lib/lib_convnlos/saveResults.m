function saveResults(rho, time, NLOSDATA, kernel)

DIR_RESULT = sprintf('results/%s', NLOSDATA.obj_name)

mkdir(DIR_RESULT);
mkdir([DIR_RESULT '/albedos']);
mkdir([DIR_RESULT '/albedos_gray']);
mkdir([DIR_RESULT '/albedos_cmap']);

methods = fieldnames(rho);

for i = 1:length(methods)
    rho = setfield(rho, methods{i}, gather(getfield(rho,  methods{i})));
end

singular_values = kernel.singular_values;
kernel_basis = kernel.basis;
result_mat = sprintf('%s/results.mat', DIR_RESULT);
NLOSDATA_save = rmfield(NLOSDATA, 'transient');
save(result_mat, 'rho', 'time', 'singular_values', 'NLOSDATA_save', 'kernel_basis');

for i = 1:length(methods)
    rho_temp = gather(getfield(rho,  methods{i}));
    saveAlbedo(rho_temp, methods{i}, [DIR_RESULT '/albedos']);
    saveAlbedo_gray(rho_temp, methods{i}, [DIR_RESULT '/albedos_gray']);
    saveAlbedo_cmap(rho_temp, methods{i}, [DIR_RESULT '/albedos_cmap']);
end



