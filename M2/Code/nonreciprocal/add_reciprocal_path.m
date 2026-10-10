function add_reciprocal_path()

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    addpath(fullfile(here, '..', 'reciprocal'));
end
