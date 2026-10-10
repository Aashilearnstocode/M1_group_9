function add_reciprocal_path()
% Puts this folder (M2/Code/nonreciprocal) and the sibling M2/Code/reciprocal folder
% (config, generate_channels, calculate_sinr, run_cga_optimizer, ...) on the MATLAB path.
% When Member 4's common framework lands, replace the second addpath with the shared folder.
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    addpath(fullfile(here, '..', 'reciprocal'));
end
