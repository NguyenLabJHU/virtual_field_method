function [s_model,s_edata] = simpleModel(model,edata, mat_surface_traction ,rParts_matrix)
% simpleModel - Updates the model struct by removing specified material blocks and their elements.
%
% model:         Original model struct.
% rParts:          Parts to be removed.
% mat_surface:   Not used in this code, but might be for surface updates.
%
% Returns:
%   s_model:      Simplified model struct with updated blocks, materials, and elements.

% Copy the full model struct to s_model
s_model = model;
s_model.all_elements = model.elements;

% List block names, their material labels, and their element ID lists
list_blockname = model.blockname;
list_blockmat = model.blockmat;
list_blockmatid = model.blockmatid;

% Identify which blocks need to be removed
rParts_flat = [rParts_matrix{:}];
idx_remove = ismember(list_blockname, rParts_flat);

% Remove the indicated blocks from all block-related lists
list_blockmat(idx_remove)    = [];
list_blockname(idx_remove)   = [];
list_blockmatid(idx_remove)  = [];

% Update s_model with the filtered block lists
s_model.blockname = list_blockname;
s_model.blockmat = list_blockmat;
s_model.blockmatid = list_blockmatid;

% --- Remove the elements and materials associated with the removed blocks ---

% Step 3: Gather matIDs that belong to the removed materials
for i_surf = 1:numel(mat_surface_traction)
    
    idxRemoveRows = [];
    matRemove = [];
    rParts = rParts_matrix{i_surf};
    
    for i = 1:numel(rParts)
        
        PartName = rParts{i};
        
        % Find the indices of blocks with this material id
        idx_Part = find(model.blockname == PartName);
    
        mat = model.blockmatid(idx_Part);
        
        % Step 4: Find rows to remove in s_model.elemmat and s_model.element
        elements_per_part = s_model.blockelem{idx_Part}(:,1);
        idx_elem = ismember(s_model.elements(:,1), elements_per_part);
        idxRemoveRows(:,i) = idx_elem;
        matRemove = [matRemove;mat];

    end

    % Total row counting to remove for this surface
    total_idx = sum(idxRemoveRows, 2);
    % Creating a new field for the removed rows
    s_model.elemmat_removed{i_surf} = s_model.elemmat(logical(total_idx),:);
    s_model.elements_removed{i_surf} = s_model.elements(logical(total_idx), :);
    s_model.matid_removed{i_surf} = matRemove;
  
    % Remove the corresponding rows
    s_model.elemmat(logical(total_idx)) = [];
    s_model.elements(logical(total_idx), :) = [];
  
    
    % Remove surface from the element
    % Get the list of element IDs to be excluded (Column 1)
    removed_ids = s_model.elements_removed{i_surf}(:, 1); 
    
    for isurf = 1:length(s_model.surfacesp)
        % Extract the current Nx6 surface matrix
        current_surf = s_model.surfacesp{isurf};
        
        % Identify rows where the 6th column (associated element ID) 
        % exists in the removed_ids list
        is_removed_element = ismember(current_surf(:, 6), removed_ids);

        % Saves the removed surfaces    
        s_model.surfacesp_removed{isurf} = current_surf(is_removed_element, :);
        
        % Filter the matrix: keep only rows where is_removed_element is FALSE
        % This ensures the surface definition matches the current mesh state
        s_model.surfacesp{isurf} = current_surf(~is_removed_element, :);
    end
        
    
end


% --- Initialization ---
surface_names = {model.all_surfaces.name};
surface_nameS_str = string(surface_names);
surfaces_to_search_str = string(mat_surface_traction);

surface_traction = struct();

% --- Main Surface Loop ---
for i = 1:numel(surfaces_to_search_str)
    surface_name_str = surfaces_to_search_str(i);
    surface_index = find(surface_nameS_str == surface_name_str);

    if isempty(surface_index)
        error(['Surface ' char(surface_name_str) ' not found in model.all_surfaces.']);
    end

    surface_data = model.all_surfaces(surface_index);
    surface_connectivity = surface_data.conn; 
    surface_traction(i).conn = surface_connectivity;

    % --- Loop through each face on the surface ---
    for iface = 1:length(surface_connectivity)
        nodal_face = surface_connectivity(iface, 2:5); 

        % Search for the parent element
        for k = 1:length(removed_ids)
            elem_nodes_temp = s_model.elements_removed{i}(k, 2:9); 
            if all(ismember(nodal_face, elem_nodes_temp))
                elem_removed_mat = s_model.elements_removed{i}(k, 1);
                surface_traction(i).element{iface} = elem_removed_mat;
                break;
            end
        end

     end
end

edata.surface_traction = surface_traction;
s_edata = edata;

end





