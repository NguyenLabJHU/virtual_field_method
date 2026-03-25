function edata_exp = accumulate_Fpre_from_edata(mydir, mymodel, gauss_order, prestress_time, matparam_complete, edata_exp, model)
% Accumulate cumulative deformation gradient (Fpre) and store in edata.steps{step}.results.eFpre_calc.
%
% Between t=0 and the user-provided prestress_time, Fpre is accumulated (product of F).
% After prestress_time, the Fpre is frozen (copied) from the prestress_time step.
%
% This version is parallel-safe: filenames are tagged with file_id to avoid collisions in parallel runs.
%
% INPUTS:
%   model         -- FE model structure (includes elements, nodes, etc)
%   edata         -- struct from preliminary_reading2 (includes steps and times)
%   gauss_order   -- integer; order for numerical integration/quadrature
%   prestress_time-- time point at which prestress is calculated and then frozen onwards
%   file_id       -- unique integer for this run (e.g. parfor loop index)
%
% OUTPUT:
%   edata         -- same struct, but with .steps{step}.results.eFpre_calc field added:
%                    Matrix [nElem, nGauss, 3, 3]: cumulative Fpre for each step

file_id = 1;

% Read the text lines from the template .feb model file as before
dataFile = fullfile(mydir,mymodel);
fid = fopen(dataFile,'r');
data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);

lines = data{1};

% Updating the material block line
lines = update_material_block_lines(lines, matparam_complete, model);

% --- Save modified lines back to a unique .feb file for this run (parallel-safe change) ---
unique_model_name = sprintf('modeltilde.feb', file_id); % <-- unique filename
full_unique_model_name = fullfile(mydir, unique_model_name);
fid = fopen(full_unique_model_name,'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);

% --- Run febio simulation with new material properties, using unique output filename ---
unique_output_name = sprintf('output.txt', file_id); % <-- unique filename
full_unique_output_name = fullfile(mydir, unique_output_name);
fprintf('FEBio running to calculate Fpre of this run (file_id=%d)\n', file_id);

t_pre = tic;
febio_bin = get_febio_path();
febio_cmd = sprintf('"%s" -i "%s" > "%s" 2>&1', ...
    febio_bin, full_unique_model_name, full_unique_output_name);

[status, cmdout] = system(febio_cmd);
elapsed = toc(t_pre);

fprintf('Computing the cumulative Fpre took %.3f seconds\n', elapsed);

%% Check for convergence
% Check for failure conditions
if status ~= 0 || contains(cmdout, 'E R R O R   T E R M I N A T I O N')
    error_flag = true;
    fprintf('FEBio terminated with error after completing timing info');
    % leave the function after the error
    return
else
    error_flag = false;
end

% --- Handle unique log file for this run ---
unique_log_name = sprintf('modeltilde.log', file_id); % <-- unique filename
full_unique_log_name = fullfile(mydir, unique_log_name);

% Model and log files for further parsing/results
mymodel = unique_model_name;
myexpdata = unique_log_name;

t_pre_read = tic;

% Parsing data from the unique output files just generated

[model, edata] = preliminary_reading(mydir, mymodel, myexpdata, edata_exp, model,prestress_time);

elapsed = toc(t_pre_read);

fprintf('Reading the cumulative Fpre took %.3f seconds\n', elapsed);

coords = model.nodes(:,2:4);                 % Get reference coordinates from model
nSteps = numel(edata.steps);                 % Total number of time steps

% ---- Find the corresponding index for the prestress time ----
% Searches for the last edata.times <= prestress_time (handles floating point times robustly)
prestress_step = find(edata.times < prestress_time, 1, 'last'); % find last time less than or equal to target
if isempty(prestress_step)
    error('Prestress time is earlier than the provided dataset.');
end

nElem = size(model.elements,1);              % Number of elements
[gauss_points, ~] = get_gauss_points(3, gauss_order); % Get Gauss integration points
nGauss = size(gauss_points,1);              % Number of gauss points per element

% Initialize "last Fpre" as identity for every element and gauss point
Fpre_last = zeros(nElem, nGauss, 3, 3);
for e = 1:nElem
    for g = 1:nGauss
        Fpre_last(e,g,:,:) = eye(3);
    end
end

% ---- Step through every time step to compute Fpre ----
for step = 1:nSteps
    if step <= prestress_step
        % --- Before (and at) prestress step: accumulate Fpre ---
        disp = edata.steps{step}.results.edisp;               % Get node displacements for this step
        
        % Compute deformation gradient at each element, gauss point (returns nElem x nGauss x 3 x 3 array)
        [Fcurr, ~] = compute_deformation_gradient(model, coords, disp, gauss_order);

        % Accumulate total Fpre by multiplying previous Fpre * current F (matrix product)
        Fpre = zeros(nElem, nGauss, 3, 3);
        for e = 1:nElem
            for g = 1:nGauss
                Fc = squeeze(Fcurr(e,g,:,:));    % Current deformation gradient at this elem/gp
                Fp = squeeze(Fpre_last(e,g,:,:));% Cumulative Fpre up to previous step
                Fpre(e,g,:,:) = Fp * Fc;         % Multiply: prev cumulative * this increment
            end
        end

        % Store cumulative Fpre inside edata struct for this step
        edata_exp.steps{step}.results.eFpre_calc = Fpre;
        
        % Save current cumulative Fpre for next step's initialization
        Fpre_last = Fpre;

        if step == prestress_step
            % Save the frozen state of Fpre for future steps
            Fpre_frozen = Fpre;
        end
        
    else
        % --- After prestress: Fpre is frozen ---
        % Store the prestress Fpre value at subsequent timesteps
        edata_exp.steps{step}.results.eFpre_calc = Fpre_frozen;
    end

    % Passing to edata_exp the information of the current Fpre simualtion    
    
    % Copy ecoords to the new ecoords_prestress field
    edata_exp.steps{step}.results.ecoords_prestress = edata.steps{step}.results.ecoords;
    
    % Copy edisp to the new edisp_prestress field
    edata_exp.steps{step}.results.edisp_prestress = edata.steps{step}.results.edisp;
end

% --- Clean up temporary files for this run ---

full_unique_xplt_name = strrep(full_unique_model_name, '.feb', '.xplt');

delete(fullfile(full_unique_model_name));      % Remove unique .feb file
delete(fullfile(full_unique_log_name));        % Remove unique .log file
delete(fullfile(full_unique_output_name));     % Remove output .txt file
delete(fullfile(full_unique_xplt_name));       % Remove output .xplt file

end