function [datau, dataF] = simulate_febio_uniform(mydir, mymodel, matparam_sweep, ...
    ground_truth_mat, nnd, nel, param_ind, changing_matrix, ops_matrix_struct, model)
%SIMULATE_FEBIO_UNIFORM Runs a single parameterized FEBio simulation.
%
%   [datau, dataF] = simulate_febio_uniform(mydir, mymodel, matparam_sweep, 
%       ground_truth_mat, nnd, nel, param_ind, changing_matrix, ops_matrix_struct, model)
%
%   Inputs:
%       mydir             - Directory containing the model template
%       mymodel           - Model template filename
%       matparam_sweep    - Matrix of parameter values for this sweep
%       ground_truth_mat  - Matrix of original (reference) material parameters
%       nnd               - Number of nodes (for storage allocation)
%       nel               - Number of elements (for storage allocation)
%       param_ind         - Index of parameter combination to run
%       changing_matrix   - Cell array indicating properties/materials to change
%       ops_matrix_struct - Structure specifying special parameter operations
%       model             - Model struct defining types of material models
%
%   Outputs:
%       datau             - Matrix of nodal displacements (nnd x 7)
%       dataF             - Matrix of element deformation/gradient data (nel x 19)
%
%   The function modifies the material block for a specific parameter combination,
%   runs a FEBio simulation, parses the results, and handles errors gracefully.

global last_time

% Preallocate storage for results
datau = zeros(nnd, 7);          % Nodal displacements
dataF = zeros(nel, 19);         % Deformation and gradient data

% Number of materials that will be modified
nummat = size(changing_matrix, 2);

% Load the model template file
full_unique_model_name = fullfile(mydir, mymodel);
fid = fopen(full_unique_model_name, 'r');
data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);
data = data{1};                 % Flatten cell array
Nline = numel(data);

% Prepare matrix of parameter values corresponding to this simulation
matparam_line = matparam_sweep(param_ind,:);
matparam = ground_truth_mat;
param2change = changing_matrix{1, param_ind};
matparam(param2change, :) = matparam_line;

% Apply any defined parameter operations before simulation
totalNparam = length(matparam(1, :));  % Number of properties per material
totalNmat = length(matparam(:, 1));    % Number of materials

for idx_Mat = 1:totalNmat
    tgt_row = ops_matrix_struct(idx_Mat).tgt_row;
    src_row = ops_matrix_struct(idx_Mat).src_row;
    A_matrix = ops_matrix_struct(idx_Mat).A;
    B_matrix = ops_matrix_struct(idx_Mat).B;
    
    src_mult = zeros(1, totalNparam);

    if isempty(A_matrix)
        continue % No special op for this material
    end

    for n_src = 1:length(src_row)
        src_mult = src_mult + matparam(src_row(n_src), :) * A_matrix{n_src};
    end

    matparam(tgt_row, :) = (src_mult + matparam(tgt_row, :) * B_matrix)';
end

% Update the material property lines within the model data
data = update_material_block_lines(data, matparam, model);

% Write a new .feb (FEBio input) file for this simulation
unique_variation_name = sprintf('modeltilde.feb'); % Generated filename
full_unique_variation_name = fullfile(mydir, unique_variation_name);
fid = fopen(full_unique_variation_name, 'w');
for I = 1:Nline
    fprintf(fid, '%s\n', char(data{I}));
end
fclose(fid);

% Call FEBio to run simulation via command line
febio_bin = get_febio_path();
unique_output_name = sprintf('output.txt'); % Generated output log
full_unique_output_name = fullfile(mydir, unique_output_name);
febio_cmd = sprintf('"%s" -i "%s" > "%s" 2>&1', ...
    febio_bin, full_unique_variation_name, full_unique_output_name);
[status, cmdout] = system(febio_cmd);

%% Check for simulation errors and handle gracefully
if status ~= 0 || contains(cmdout, 'E R R O R   T E R M I N A T I O N')
    error_flag = true;
    disp('FEBio terminated with error after completing timing info');
    datau = nan(nnd, 7);   % Fill results with NaN on error
    dataF = nan(nel, 19);
    return % Exit on error
else
    error_flag = false;
end

%% Parse the output .log file from FEBio to extract results
% Open output log file for reading
unique_log_name = sprintf('modeltilde.log'); % The output log file name
full_unique_log_name = fullfile(mydir, unique_log_name);
fres = fopen(full_unique_log_name, 'r');
if fres == -1
    error('Could not open log file.');
end

% Search for the time step header corresponding to last_time
time_string = sprintf('Time = %.6g', last_time);
s = fgetl(fres);
found_time = false;
while ischar(s)
    % Robust search: matches times formatted as floats (e.g., 1.5 == 1.500)
    if contains(s, 'Time =')
        tmp_time = sscanf(s, 'Time = %f');
        if ~isempty(tmp_time) && abs(tmp_time - last_time) < 1e-8
            found_time = true;
            break;
        end
    end
    s = fgetl(fres);
end

if ~found_time
    fclose(fres);
    error('Requested time %g not found in log file.', last_time);
end

% Search for Displacement Data header line
displ_header = 'Data = x;y;z;ux;uy;uz';
while ischar(s)
    if strcmp(s, displ_header)
        break;
    end
    s = fgetl(fres);
end
if ~ischar(s)
    fclose(fres);
    error('Displacement data header not found.');
end

% Parse the displacement data block
datau = zeros(nnd, 7);
for idk = 1:nnd
    vals = fscanf(fres,'%i %f %f %f %f %f %f\n',[1,7]);
    if isempty(vals) || numel(vals) ~= 7
        fclose(fres);
        error('Displacement data incomplete at node %d.', idk);
    end
    datau(idk, :) = vals;
end

fclose(fres);  % Clean up log file handle

end