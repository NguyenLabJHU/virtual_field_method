function [ground_truth_mat, matparam_sweep, matparam_complete] = SweepMatrix(model, changing_matrix, x, Normalizer, ops_matrix_struct)
% SweepMatrix creates matrices to be used during the generation of virtual fields (matparam_sweep),
% groun_truth_mat (matrix of ground truth parameter), matparam_complete that is used to generate 
% the prestress data from the guess in the optmization routine.
%
%
%   Inputs:
%       model              - Struct containing material models and properties
%       changing_matrix    - Cell array describing which properties to change in which material
%       x                  - Vector of input parameters to be swept
%       Normalizer         - Normalization factor for the parameters
%       ops_matrix_struct  - Struct describing how to apply parameter combinations/operations
%
%   Outputs:
%       ground_truth_mat   - Unmodified material parameter matrix (original)
%       matparam_sweep     - Matrix showing parameter sweep modifications
%       matparam_complete  - Matrix showing all modifications and applied operations
%
%   The function applies new values to selected material properties and manages 
%   relationships/constraints between properties (e.g., incompressibility).

    ground_truth_mat = model.matprop;   % Store the original material property matrix
    matparam_complete = ground_truth_mat;  % Will hold final modified material parameters
    corresponding = cell2mat(changing_matrix(1,:));   % Indices of materials to change
    parameters = x .* Normalizer;       % Normalize parameters to correct scale
    matparam_sweep = zeros([length(corresponding), size(ground_truth_mat, 2)]);

    % Initialize material parameter matrix to be updated
    for k = 1:max(corresponding)
        % Only update materials listed in corresponding
        if ismember(k, corresponding)
            idx_int = find(corresponding == k, 1, 'first');
            idx_last = find(corresponding == k, 1, 'last');
            nParam = idx_last - idx_int + 1;

            if idx_int - idx_last == 0  % Single parameter for this material
                % Special case for retina (material index 3)
                if k == 3
                    matparam(k,1) = parameters(idx_int);
                    matparam(k,2) = parameters(idx_int) * 20;
                else
                    matparam(k,1) = parameters(idx_int);
                    matparam(k,2) = parameters(idx_int) * 1000;
                end
            else   % Multiple parameters for this material
                for iParam = 1:nParam
                    matparam(k,iParam) = parameters(idx_int + iParam - 1);
                end
            end
        end
    end

    % Assign parameter values based on changing_matrix for each property change
    for col = 1:size(changing_matrix, 2)
        mat_idx = changing_matrix{1, col};                        % Material index
        mat_type = lower(strtrim(model.matmodel(mat_idx)));       % Material model type
        mat_change = changing_matrix{4, col};                     % Parameter name to change

        % Get the property index for the material type and parameter name
        prop = mat_change2prop(mat_type, mat_change);

        col_param_prop = prop;
        row_param = changing_matrix{1, col};
        col_param = changing_matrix{2, col};
        matparam_sweep(col, :) = ground_truth_mat(row_param, :);   % Start from ground truth
        matparam_sweep(col, col_param_prop) = matparam(row_param, col_param);  % Update swept param

        % Update complete matrix with new parameter value
        matparam_complete(mat_idx, col_param_prop) = matparam(row_param, col_param);
    end

    % Apply any specified linear operations to parameters as defined in ops_matrix_struct
    totalNparam = length(matparam_complete(1, :));  % Number of parameters per material
    totalNmat = length(matparam_complete(:,1));     % Number of materials

    for idx_Mat = 1:totalNmat
        tgt_row = ops_matrix_struct(idx_Mat).tgt_row;   % Target material to update
        src_row = ops_matrix_struct(idx_Mat).src_row;   % Source materials for operation
        A_matrix = ops_matrix_struct(idx_Mat).A;        % Multipliers for sources
        B_matrix = ops_matrix_struct(idx_Mat).B;        % Multiplier for target material
        
        src_mult = zeros(1, totalNparam);               % Initialize summed parameter vector

        if isempty(A_matrix)
            continue  % No operation for this material
        end

        % Sum over all sources defined for this operation
        for n_src = 1:length(src_row)
            src_mult = src_mult + matparam_complete(src_row(n_src), :) * A_matrix{n_src};
        end

        % Update the target material parameters with operation result
        matparam_complete(tgt_row, :) = (src_mult + matparam_complete(tgt_row, :) * B_matrix)';
    end

    % Enforce incompressibility: For materials without a prescribed K parameter,
    % set K = 1000 * first parameter (often c1)
    param_names = changing_matrix(4,:);      % Names of changed parameters
    materials_number = changing_matrix(1,:); % Material indices for those parameters

    % Find materials for which 'k' or 'K' is explicitly changed
    idx_k = cellfun(@(x) ~isempty(regexpi(x, '^k$', 'once')), param_names);
    materials2exclude = cell2mat(materials_number(idx_k));

    % Loop over all materials - update K if not explicitly set
    for i = 1:totalNmat
        if ismember(i, materials2exclude)
            continue % Skip materials where K was set
        else
            matparam_complete(i,3) = matparam_complete(i,1) * 1000; % Set K parameter
        end
    end

end