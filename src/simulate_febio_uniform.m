function [datau,dataF]=simulate_febio_uniform(mydir,mymodel,matparam_sweep,...
    ground_truth_mat,nnd,nel,param_ind,changing_matrix,ops_matrix_struct,model)

global last_time
 
%dimension storage for simulation results
%datau=storage for displacement (u)
datau = zeros(nnd,7);

%dataF=storage for deformation gradient (F) and discplacement gradient (G)
dataF = zeros(nel,19);

%simulate data with target properties. It reads in template.
%nummat=Number of materials to modify
nummat = size(changing_matrix,2);

% Read file as before
full_unique_model_name = fullfile(mydir, mymodel);
fid = fopen(full_unique_model_name,'r');
data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);
data = data{1}; % flatten cell array
Nline = numel(data);

matparam_line = matparam_sweep(param_ind,:);
matparam = ground_truth_mat;
param2change = changing_matrix{1,param_ind};
matparam(param2change,:) = matparam_line;

%OBS%
%THIS PART IS HARDCODED TO ADDRESS THE SCLERA BEING SPLIT IN THREE%

% Performing operations between parameters
totalNparam =length(matparam(1,:));
%Changing the K column from the matparam_complete
totalNmat = length(matparam(:,1));

for idx_Mat = 1:totalNmat
    tgt_row = ops_matrix_struct(idx_Mat).tgt_row;
    src_row = ops_matrix_struct(idx_Mat).src_row;
    A_matrix = ops_matrix_struct(idx_Mat).A;
    B_matrix = ops_matrix_struct(idx_Mat).B;

    src_mult = zeros(1,totalNparam);

    if isempty(A_matrix)
        continue
    end

    for n_src = 1:length(src_row)
    
       src_mult = src_mult + matparam(src_row(n_src),:)* A_matrix{n_src}; 
    
    end

    matparam(tgt_row,:) = (src_mult + matparam(tgt_row,:) * B_matrix)';
end





data = update_material_block_lines(data, matparam, model);

unique_variation_name = sprintf('modeltilde.feb'); % <-- unique filename
full_unique_variation_name = fullfile(mydir, unique_variation_name);

fid = fopen(full_unique_variation_name, 'w');
for I = 1:Nline
    fprintf(fid, '%s\n', char(data{I}));
end
fclose(fid);

%run febio simulations with new material properties.
%status 0=Success (FEBio ran without errors)
%status 1=Generic error (e.g., Convergence problem)
%status -1=Command could not be executed (e.g., path incorrect)
febio_bin = get_febio_path();

unique_output_name = sprintf('output.txt'); % <-- unique filename
full_unique_output_name = fullfile(mydir, unique_output_name);


febio_cmd = sprintf('"%s" -i "%s" > "%s" 2>&1', ...
    febio_bin, full_unique_variation_name, full_unique_output_name);

[status, cmdout] = system(febio_cmd);

%% Check for convergence
% Check for failure conditions
if status ~= 0 || contains(cmdout, 'E R R O R   T E R M I N A T I O N')
    error_flag = true;
    disp('FEBio terminated with error after completing timing info');
    
    % Match original nodedat dimensions but fill with NaN
    datau = nan(nnd, 7);    
    
    % Match original elemdat dimensions but fill with NaN
    dataF = nan(nel, 19);   
    
    %leave the function after the error
    return
else
    error_flag = false;
end


%% Parse FEBio's log file o extract nodal displacement and deformation gradients

% Open file for reading
unique_log_name = sprintf('modeltilde.log'); % <-- unique filename
full_unique_log_name = fullfile(mydir, unique_log_name);

fres = fopen(full_unique_log_name,'r');
if fres == -1
    error('Could not open log file.');
end

% Construct string to search for Time
time_string = sprintf('Time = %.6g', last_time);

% -- Search for the Time header --
s = fgetl(fres);
found_time = false;
while ischar(s)
    % Robust search: also matches slightly formatted times (e.g., 1.5 == 1.500)
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

% -- Search for Displacement Data header --
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

% -- Parse displacement data --
datau = zeros(nnd, 7);
for idk = 1:nnd
    vals = fscanf(fres,'%i %f %f %f %f %f %f\n',[1,7]);
    if isempty(vals) || numel(vals) ~= 7
        fclose(fres);
        error('Displacement data incomplete at node %d.', idk);
    end
    datau(idk, :) = vals;
end


fclose(fres);

