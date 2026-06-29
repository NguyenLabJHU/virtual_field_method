function [IVW, EVW, TieVW, cost_func] = calc_virtual_work_variation_integration_simp( ...
    path, mymodel, model, edata2, matparam, matparam_sweep, ground_truth_mat, ...
    p_app, gauss_order, eps, changing_matrix, ops_matrix_struct, full_model, nodal_force)

%--------------------------------------------------------------------------
% calc_virtual_work_variation_integration_simp
%--------------------------------------------------------------------------
% Computes:
%   IVW       = Internal Virtual Work
%   EVW       = External Virtual Work from pressure/surface loads
%   TieVW     = Virtual Work contribution from tie-contact penalty forces
%   nodal_VW  = Virtual Work contribution from prescribed nodal forces
%   cost_func = normalized residual between internal and external virtual work
%
% This function is used inside material-parameter optimization. For each
% parameter set, it either:
%   1) runs FEBio to compute the virtual displacement field, or
%   2) loads cached nodal/element outputs from CSV files.
%
% Then it computes the virtual strain field, integrates stress : virtual strain
% over the volume, computes external virtual work contributions, and accumulates
% the normalized cost.
%
% INPUTS:
%   path
%       Struct containing folders. Expected fields:
%           path.data : folder containing FEBio model/data
%           path.VF   : folder used for Virtual Field cache files
%
%   mymodel
%       FEBio model filename, for example 'model.feb'.
%
%   model
%       Current/reduced model struct. Contains nodes, elements, material IDs,
%       pressure surfaces, tie-contact information, etc.
%
%   edata2
%       FEBio/experimental result struct from preprocessing. Contains time
%       steps, nodal coordinates, displacements, deformation gradients,
%       prestress quantities, and tie-contact integration-point maps.
%
%   matparam
%       Matrix of material parameters indexed by material number.
%       Example: matparam(matnum,:) gives parameters for material matnum.
%
%   matparam_sweep
%       Material parameter set/sweep passed to FEBio simulation.
%
%   ground_truth_mat
%       Ground-truth/reference material parameters used by simulation routine.
%
%   p_app
%       Applied pressure/load values for each pressure surface.
%
%   gauss_order
%       Quadrature order for volume and surface integration.
%
%   eps
%       Penalty stiffness/parameter used in tie-contact virtual work.
%
%   changing_matrix
%       Matrix defining which material parameters are being varied.
%       Number of columns = number of parameter sets, nParam.
%
%   ops_matrix_struct
%       Auxiliary structure used by simulate_febio_uniform.
%
%   full_model
%       Full model struct. Used because some arrays, such as F_all_pre and
%       delta_e, may be indexed according to the full model element ordering.
%
%   nodal_force
%       Nodal force array used for nodal virtual work calculation.
%
% OUTPUTS:
%   IVW
%       Internal Virtual Work for the last parameter set evaluated.
%
%   EVW
%       Pressure/surface External Virtual Work for the last parameter set.
%
%   TieVW
%       Tie-contact Virtual Work for the last parameter set.
%
%   cost_func
%       Final cost value accumulated over all parameter sets:
%
%           sqrt(sum(((IVW - EVW - TieVW - nodal_VW)/EVW_denom)^2))
%
% GLOBALS:
%   totalRunCount
%       Indicates whether this is the first run. If first run, FEBio is run
%       and caches are written. Otherwise, cached FEBio outputs are loaded.
%
%   ForwardCount
%       Used to decide whether EVW/TieVW/Nodal caches should be reused.
%
%   EVW_vector
%       Stores individual pressure/surface EVW contributions.
%
%   TieVW_vector
%       Stores individual tie-contact face contributions.
%
%   prestress_time
%       Time defining whether prestress/Fpre should be included.
%--------------------------------------------------------------------------

global totalRunCount ForwardCount EVW_vector TieVW_vector prestress_time


% Kick off parallel pool
if isempty(gcp('nocreate'))

    N = str2double(getenv('SLURM_CPUS_PER_TASK'));
    
    if isnan(N) || N == 0
        N = feature('numcores'); % fallback se rodar local
    end
    
    parpool('local', 8);

end


%% ------------------------------------------------------------------------
%  Basic mesh and time-step setup
% -------------------------------------------------------------------------

% Number of volume elements and nodes in the current/reduced model
nel = length(model.elements);
nnd = length(model.nodes);

% Time-step information from FEBio/experimental data
times           = edata2.times;
first_time_step = 1;              % Reference configuration
last_time_step  = length(times);  % Final/deformed configuration

% Time step used as the beginning of the elastic/prestress phase
prestress_step = find(edata2.times <= prestress_time, 1, 'last');

% Reference configuration at the first time step
% ecoords0: nodal coordinates in the reference configuration
% edisp0  : nodal displacements at reference time
res0     = edata2.steps{1, first_time_step}.results;
ecoords0 = res0.ecoords;
edisp0   = res0.edisp;

% Final/deformed configuration
% ecoords: nodal coordinates at final time
% edisp  : nodal displacements at final time
resN    = edata2.steps{1, last_time_step}.results;
ecoords = resN.ecoords;
edisp   = resN.edisp;

% Previous time step, used mainly for tie-contact/contact calculations
resN_1    = edata2.steps{1, last_time_step - 1}.results;
ecoords_1 = resN_1.ecoords;
edisp_1   = resN_1.edisp;

% Configuration at the beginning of the elastic/prestress phase
resN_pre0    = edata2.steps{1, prestress_step}.results;
ecoords_pre0 = resN_pre0.ecoords;
edisp_pre0   = resN_pre0.edisp;

% Prestress-related configuration at the final time
% ecoords_pre: coordinates after applying prestress
% edisp_pre  : prestress displacement used as reference for virtual field
resN_pre    = edata2.steps{1, last_time_step}.results;
ecoords_pre = resN_pre.ecoords_prestress;
edisp_pre   = resN_pre.edisp_prestress;

% Gauss order used for hexahedral volume integration
gauss_order_hex = gauss_order;

% Model name without .feb extension, used for cache filenames
modelName = erase(mymodel, '.feb');


%% ------------------------------------------------------------------------
%  Cache file paths
% -------------------------------------------------------------------------

% Cache files for FEBio virtual-field outputs.
% nodedat cache contains nodal output for all parameter sets.
% elemdat cache contains element output for all parameter sets.
nodeCachePath = fullfile(path.VF, ['nodedat_cached_', modelName, '.csv']);
elemCachePath = fullfile(path.VF, ['elemdat_cached_', modelName, '.csv']);

% Cache files for virtual work contributions.
% energyFile stores pressure/surface EVW.
% tieFile stores tie-contact virtual work.
% nodalFile stores nodal-force virtual work.
energyFile = fullfile(path.VF, ['EVW_cached_simp_', modelName, '.csv']);
tieFile    = fullfile(path.VF, ['TieVW_cached_simp_', modelName, '.csv']);
nodalFile  = fullfile(path.VF, ['Nodal_VW_cached_simp_', modelName, '.csv']);

% Original EVW cache, used as cost-function normalization denominator
% when available.
EVW_original_file = fullfile(path.VF, ['EVW_cached_', modelName, '.csv']);

% Cache file for virtual strain tensor delta_e for each parameter set.
% Stored as:
%   delta_e_struct.param_1
%   delta_e_struct.param_2
%   ...
cachefile_delta_e = ['delta_e_struct_', modelName, '.mat'];
fulldelta_e_File  = fullfile(path.VF, cachefile_delta_e);


%% ------------------------------------------------------------------------
%  Precompute deformation gradient and full-model element map
% -------------------------------------------------------------------------

% F_all contains deformation gradient Fbar for the current/reduced model.
% Expected size:
%   [nel x nGauss x 3 x 3]
[F_all, ~] = compute_deformation_gradient(model, ecoords0, edisp, gauss_order_hex);

% F_all_pre contains the prestress/predeformation gradient from the full model.
% It is only used if prestress_time is nonzero.
% This array may be indexed using full_model element ordering.
if abs(prestress_time) > 1e-6
    F_all_pre = edata2.steps{1, last_time_step}.results.eFpre_calc;
else
    F_all_pre = [];
end

% Map each element in model.elements to its row in full_model.elements.
% For current/reduced model element k:
%
%   kk_full = map_to_full(k)
%
% gives the corresponding row index in full_model.
% This is needed because arrays such as F_all_pre and delta_e may follow
% full_model ordering.
[~, map_to_full] = ismember(model.elements(:, 1), full_model.elements(:, 1));

if any(map_to_full == 0)
    error('Some model.elements were not found in full_model.elements');
end


%% ------------------------------------------------------------------------
%  Parameter and output initialization
% -------------------------------------------------------------------------

% Number of material-parameter sets to evaluate
nParam = size(changing_matrix, 2);

% Accumulated cost over all parameter sets
cost_func = 0;

% Output arrays for caching FEBio virtual-field data.
% nodedat has 7 columns per parameter set:
%   [node_id, x, y, z, ux, uy, uz]
%
% elemdat has 19 columns per parameter set:
%   [element_id, ... deformation/prestress quantities ...]
nodedat_out = zeros(nnd, 7  * nParam);
elemdat_out = zeros(nel, 19 * nParam);

% Output arrays for caching scalar virtual work values for each parameter set
EVW_out      = zeros(nParam, 1);
TieVW_out    = zeros(nParam, 1);
nodal_VW_out = zeros(nParam, 1);

% Initialize outputs in case the function exits early
IVW   = 0;
EVW   = 0;
TieVW = 0;


%% ------------------------------------------------------------------------
%  Load or compute FEBio virtual-field simulation data
% -------------------------------------------------------------------------

if totalRunCount == 1

    % First run: run FEBio simulations and later write cache files
    NewVirtualWorkFlag = 1;

else

    % Later runs: load previously computed FEBio outputs from cache
    tic
    nodedat_cached = readmatrix(nodeCachePath);
    elemdat_cached = readmatrix(elemCachePath);
    fprintf('reading nodedat_cached/elemdat_cached took %.4f s\n', toc);

    NewVirtualWorkFlag = 2;

end


%% ------------------------------------------------------------------------
%  Load EVW/TieVW/Nodal caches if available
% -------------------------------------------------------------------------

% Cache flags. If ForwardCount == 1, force recomputation/overwrite.
hasEnergyCache = isfile(energyFile) && ForwardCount ~= 1;
hasTieCache    = isfile(tieFile)    && ForwardCount ~= 1;
hasNodalCache  = isfile(nodalFile)  && ForwardCount ~= 1;

if hasEnergyCache
    EVW_out = readmatrix(energyFile);
end

if hasTieCache
    TieVW_out = readmatrix(tieFile);
end

if hasNodalCache
    nodal_VW_out = readmatrix(nodalFile);
end

% Original EVW used as denominator in the cost function if available
if isfile(EVW_original_file)
    EVW_original = readmatrix(EVW_original_file);
    hasEVWOriginal = true;
else
    EVW_original = [];
    hasEVWOriginal = false;
end


%% ------------------------------------------------------------------------
%  Gauss points
% -------------------------------------------------------------------------

% Volume integration: 3D hexahedral elements
dimension_hex = 3;
[gauss_points_elem, weights_elem] = get_gauss_points(dimension_hex, gauss_order_hex);

% Surface integration: 2D quadrilateral faces
gauss_order_surface = gauss_order;
dimension_surface   = 2;
[gauss_points_surf, weights_surf] = get_gauss_points(dimension_surface, gauss_order_surface);


%% ------------------------------------------------------------------------
%  Parameter sweep
% -------------------------------------------------------------------------

for param_ind = 1:nParam

    fprintf('\n--- Parameter %d/%d ---\n', param_ind, nParam);


    %% --------------------------------------------------------------------
    %  Load or simulate virtual displacement field
    % ---------------------------------------------------------------------

    if NewVirtualWorkFlag == 1

        % First run: run FEBio for this parameter set
        mydir = path.data;

        tic
        [nodedat, elemdat] = simulate_febio_uniform( ...
            mydir, mymodel, matparam_sweep, ground_truth_mat, nnd, nel, ...
            param_ind, changing_matrix, ops_matrix_struct, model);
        fprintf('simulate_febio_uniform took %.4f s\n', toc);

        % Virtual displacement = simulated displacement - prestress reference
        nodedat(:, 5:7) = nodedat(:, 5:7) - edisp_pre;

        % Store this parameter set in output cache arrays
        nodedat_out(:, (param_ind - 1) * 7  + 1 : param_ind * 7)  = nodedat;
        elemdat_out(:, (param_ind - 1) * 19 + 1 : param_ind * 19) = elemdat;

    else

        % Later runs: extract this parameter set from cached matrices
        nodedat = nodedat_cached(:, (param_ind - 1) * 7  + 1 : param_ind * 7);
        elemdat = elemdat_cached(:, (param_ind - 1) * 19 + 1 : param_ind * 19);

    end


    %% --------------------------------------------------------------------
    %  NaN check
    % ---------------------------------------------------------------------

    % If FEBio failed or output is corrupted, abort this objective evaluation
    if any(isnan(nodedat(:))) || any(isnan(elemdat(:)))
        IVW       = NaN;
        EVW       = NaN;
        TieVW     = NaN;
        cost_func = NaN;
        return;
    end


    %% --------------------------------------------------------------------
    %  Virtual displacement and virtual strain
    % ---------------------------------------------------------------------

    tic

    % Virtual displacement field at nodes
    % nodedat columns 5:7 are ux, uy, uz after subtracting edisp_pre
    delu = nodedat(:, 5:7);

    % Field name used inside delta_e_struct cache
    fieldname = sprintf('param_%d', param_ind);

    need_to_compute = false;
    delta_e_struct  = struct();

    % Load virtual strain tensor delta_e from cache if available.
    % Otherwise compute it and save it.
    if exist(fulldelta_e_File, 'file')

        S = load(fulldelta_e_File);

        if isfield(S, 'delta_e_struct') && isfield(S.delta_e_struct, fieldname)

            % Virtual strain tensor for this parameter set
            % Expected size:
            %   [nElem_full_or_model x nGauss x 3 x 3]
            delta_e = S.delta_e_struct.(fieldname);

        else

            need_to_compute = true;

            if isfield(S, 'delta_e_struct')
                delta_e_struct = S.delta_e_struct;
            else
                delta_e_struct = struct();
            end

        end

    else

        need_to_compute = true;
        delta_e_struct  = struct();

    end

    if need_to_compute

        % Compute virtual strain field at Gauss points
        [delta_e, ~] = compute_virtual_strain_integration( ...
            model, ecoords, delu, gauss_order_hex);

        fprintf('Computed new delta_e for %s\n', fieldname);

        % Save/update delta_e cache
        delta_e_struct.(fieldname) = delta_e;
        save(fulldelta_e_File, 'delta_e_struct');

    end


    %% --------------------------------------------------------------------
    %  Internal Virtual Work
    % ---------------------------------------------------------------------
    % IVW = integral_over_volume sigma : delta_e dV
    %
    % element_IVW_array stores each element contribution. The sum after the
    % parfor gives total IVW.

    element_IVW_array = zeros(nel, 1);

    parfor k = 1:nel

        % Element ID in current/reduced model
        ielem = model.elements(k, 1);

        % Material number assigned to this element
        matnum = model.elemmat(k);

        % Corresponding row in full_model. Needed for arrays indexed by
        % full_model ordering, such as F_all_pre and delta_e.
        kk_full = map_to_full(k);

        % Material parameters for this element/material
        imatprop = zeros(size(model.matprop, 2), 1);
        imatprop(:) = matparam(matnum, :);

        % Coordinates of the 8 nodes of the hexahedral element
        X = ecoords(model.elements(k, 2:9), :);

        % Element-level accumulators
        element_IVW        = 0;
        element_Vol        = 0;
        element_vf_norm_sq = 0;

        % Loop over Gauss points in the hexahedral element
        for gp = 1:size(gauss_points_elem, 1)

            xi     = gauss_points_elem(gp, 1);
            eta    = gauss_points_elem(gp, 2);
            zeta   = gauss_points_elem(gp, 3);
            weight = weights_elem(gp);

            % Hexahedral shape function derivatives with respect to
            % natural coordinates xi, eta, zeta.
            dN_dxi = 0.125 * [ ...
                -(1 - eta) * (1 - zeta);
                 (1 - eta) * (1 - zeta);
                 (1 + eta) * (1 - zeta);
                -(1 + eta) * (1 - zeta);
                -(1 - eta) * (1 + zeta);
                 (1 - eta) * (1 + zeta);
                 (1 + eta) * (1 + zeta);
                -(1 + eta) * (1 + zeta)];

            dN_deta = 0.125 * [ ...
                -(1 - xi) * (1 - zeta);
                -(1 + xi) * (1 - zeta);
                 (1 + xi) * (1 - zeta);
                 (1 - xi) * (1 - zeta);
                -(1 - xi) * (1 + zeta);
                -(1 + xi) * (1 + zeta);
                 (1 + xi) * (1 + zeta);
                 (1 - xi) * (1 + zeta)];

            dN_dzeta = 0.125 * [ ...
                -(1 - xi) * (1 - eta);
                -(1 + xi) * (1 - eta);
                -(1 + xi) * (1 + eta);
                -(1 - xi) * (1 + eta);
                 (1 - xi) * (1 - eta);
                 (1 + xi) * (1 - eta);
                 (1 + xi) * (1 + eta);
                 (1 - xi) * (1 + eta)];

            % Jacobian matrix from natural to physical coordinates
            dX_dxi   = dN_dxi'   * X;
            dX_deta  = dN_deta'  * X;
            dX_dzeta = dN_dzeta' * X;

            Jac_matrix = [dX_dxi; dX_deta; dX_dzeta];

            jac_det      = det(Jac_matrix);
            jac_weighted = jac_det * weight;

            % Elastic/mechanical deformation gradient from current model
            mFbar = squeeze(F_all(k, gp, :, :));

            % Prestress/predeformation gradient from full model, if used
            if abs(prestress_time) > 1e-6
                mFpre = squeeze(F_all_pre(kk_full, gp, :, :));
            else
                mFpre = eye(3);
            end

            % Total deformation gradient
            F = mFbar * mFpre;

            % Virtual strain tensor at this element/Gauss point.
            % Indexed using full-model ordering.
            dgum = squeeze(delta_e(kk_full, gp, :, :));

            % Cauchy stress tensor at this Gauss point
            sigma = computeStress(F, model, matnum, imatprop, ielem, gp);

            % sigma : delta_e = sum_ij sigma_ij * delta_e_ij
            temp = sum(sum(dgum .* sigma));

            % Add weighted Gauss-point contribution to element IVW
            element_IVW = element_IVW + temp * jac_weighted;

            % Element volume accumulator, kept for possible future use
            element_Vol = element_Vol + jac_weighted;

            % Virtual strain norm accumulator, kept for possible future use
            vf_q = sum(sum(dgum .* dgum));
            element_vf_norm_sq = element_vf_norm_sq + vf_q * jac_weighted;

        end

        % Store element contribution. Sum is done after parfor.
        element_IVW_array(k) = element_IVW;

    end

    % Total internal virtual work
    IVW = sum(element_IVW_array);


    %% --------------------------------------------------------------------
    %  External Virtual Work - pressure/surface contribution
    % ---------------------------------------------------------------------
    % EVW = integral_over_loaded_surface p * delta_u dot n dA

    if hasEnergyCache

        % Load pressure/surface EVW from cache
        EVW = EVW_out(param_ind);

    else

        EVW   = 0;
        AreaS = 0;

        % Number of pressure/load surfaces
        nLoad = length(p_app);

        for nLoad_idx = 1:nLoad

            % Number of surface elements for this load
            nels_surface = length(model.surfacesp{1, nLoad_idx});

            for ksurf = 1:nels_surface

                % Coordinates and virtual displacements of the 4 surface nodes
                X = zeros(4, 3);
                delu_nodes = zeros(4, 3);

                for node_idx = 1:4
                    node_id = model.surfacesp{1, nLoad_idx}(ksurf, node_idx + 1);

                    X(node_idx, :)          = ecoords(node_id, :);
                    delu_nodes(node_idx, :) = delu(node_id, :);
                end

                surface_EVW  = 0;
                surface_Area = 0;

                % Surface Gauss integration
                for gp = 1:size(gauss_points_surf, 1)

                    xi     = gauss_points_surf(gp, 1);
                    eta    = gauss_points_surf(gp, 2);
                    weight = weights_surf(gp);

                    % Bilinear quad shape functions
                    N = 0.25 * [ ...
                        (1 - xi) * (1 - eta);
                        (1 + xi) * (1 - eta);
                        (1 + xi) * (1 + eta);
                        (1 - xi) * (1 + eta)];

                    dN_dxi = 0.25 * [ ...
                        -(1 - eta);
                         (1 - eta);
                         (1 + eta);
                        -(1 + eta)];

                    dN_deta = 0.25 * [ ...
                        -(1 - xi);
                        -(1 + xi);
                         (1 + xi);
                         (1 - xi)];

                    % Surface Jacobian and normal
                    dX_dxi  = dN_dxi'  * X;
                    dX_deta = dN_deta' * X;

                    area_vec = cross(dX_dxi, dX_deta);
                    jac_surf = norm(area_vec);

                    jac_weighted = jac_surf * weight;
                    normal       = area_vec / jac_surf;

                    % Interpolated virtual displacement at this surface GP
                    delum = N' * delu_nodes;

                    % Pressure work contribution
                    surface_EVW = surface_EVW + ...
                        p_app(nLoad_idx) * dot(delum, normal) * jac_weighted;

                    surface_Area = surface_Area + jac_weighted;

                end

                EVW = EVW + surface_EVW;

                % Store per-surface-element contribution in global vector
                EVW_vector(nLoad_idx, ksurf, param_ind) = surface_EVW;

                AreaS = AreaS + surface_Area;

            end

        end

        EVW_out(param_ind) = EVW;

    end


    %% --------------------------------------------------------------------
    %  Tie-contact Virtual Work
    % ---------------------------------------------------------------------
    % TieVW is computed using tie-contact integration point maps and penalty
    % forces based on the contact gap difference.

    if hasTieCache

        TieVW = TieVW_out(param_ind);

    else

        TieVW  = 0;
        TieArea = 0;

        % Kept for compatibility/debugging with previous version
        test = 1;
        diff_matrix = [];
        gap_info = [];

        % Only compute if tie-contact information exists
        if isfield(model, 'tie_contact') && ~isempty(model.tie_contact)

            nfaces = size(model.tie_contact.primary_nodes, 1);

            for iface = 1:nfaces

                % Integration point map for this primary face
                ip_map0 = edata2.all_ip_map{1, 1}{iface, 1};

                % Primary face node IDs
                primary_face_nodes = model.tie_contact.primary_nodes(iface, :);

                surface_TieVW = 0;
                surface_Area  = 0;

                % Primary face coordinates and virtual displacements
                X_face_current  = zeros(4, 3);
                X_face_1        = zeros(4, 3);
                X_face_ref      = zeros(4, 3);
                delu_nodes_prim = zeros(4, 3);

                for kface = 1:4

                    idx = find(model.nodes(:, 1) == primary_face_nodes(kface + 1));

                    X_face_current(kface, :)  = ecoords_pre(idx, 1:3);
                    X_face_1(kface, :)        = ecoords_1(idx, 1:3);
                    X_face_ref(kface, :)      = ecoords0(idx, 1:3);
                    delu_nodes_prim(kface, :) = delu(primary_face_nodes(kface + 1), :);

                end

                % Loop over tie-contact integration points
                for ip = 1:length(ip_map0)

                    % Primary integration point coordinates
                    xi  = ip_map0(ip).prim_rs(1);
                    eta = ip_map0(ip).prim_rs(2);

                    % Primary quad shape functions
                    Np = 0.25 * [ ...
                        (1 - xi) * (1 - eta);
                        (1 + xi) * (1 - eta);
                        (1 + xi) * (1 + eta);
                        (1 - xi) * (1 + eta)];

                    % Primary integration point positions
                    prim_xyz_ref      = Np' * X_face_ref;
                    prim_xyz_current  = Np' * X_face_1;
                    prim_xyz_current2 = Np' * X_face_current;

                    % Primary virtual displacement at integration point
                    delu_prim_gp = Np' * delu_nodes_prim;

                    % Secondary element/face information from the IP map
                    sec_elem_id    = ip_map0(ip).sec_elem_id;
                    sec_elem_nodes = ip_map0(ip).sec_face;

                    % Secondary face coordinates and virtual displacements
                    delu_nodes_sec = zeros(4, 3);
                    X_sec_1        = zeros(4, 3);
                    X_sec_ref      = zeros(4, 3);
                    X_sec_current  = zeros(4, 3);

                    for k2 = 1:4

                        idx2 = find(model.nodes(:, 1) == sec_elem_nodes(k2));

                        delu_nodes_sec(k2, :) = delu(sec_elem_nodes(k2), :);
                        X_sec_1(k2, :)        = ecoords_1(idx2, 1:3);
                        X_sec_current(k2, :)  = ecoords_pre(idx2, 1:3);
                        X_sec_ref(k2, :)      = ecoords0(idx2, 1:3);

                    end

                    % Secondary integration point coordinates
                    xi_sec  = ip_map0(ip).sec_rs(1);
                    eta_sec = ip_map0(ip).sec_rs(2);

                    % Secondary quad shape functions
                    N_sec = 0.25 * [ ...
                        (1 - xi_sec) * (1 - eta_sec);
                        (1 + xi_sec) * (1 - eta_sec);
                        (1 + xi_sec) * (1 + eta_sec);
                        (1 - xi_sec) * (1 + eta_sec)];

                    % Secondary integration point positions
                    sec_xyz_ref      = N_sec' * X_sec_ref;
                    sec_xyz_current  = N_sec' * X_sec_1;
                    sec_xyz_current2 = N_sec' * X_sec_current;

                    % Secondary virtual displacement at integration point
                    delu_sec_gp = N_sec' * delu_nodes_sec;

                    % Virtual gap between primary and secondary virtual displacements
                    delta_g = delu_prim_gp - delu_sec_gp;

                    % Surface Jacobian for primary face
                    dN_dxi = 0.25 * [ ...
                        -(1 - eta);
                         (1 - eta);
                         (1 + eta);
                        -(1 + eta)];

                    dN_deta = 0.25 * [ ...
                        -(1 - xi);
                        -(1 + xi);
                         (1 + xi);
                         (1 - xi)];

                    dX_dxi  = dN_dxi'  * X_face_current;
                    dX_deta = dN_deta' * X_face_current;

                    area_vec = cross(dX_dxi, dX_deta);
                    jac_surf = norm(area_vec);

                    % Integration weight from tie-contact IP map
                    w_ip = ip_map0(ip).weight;
                    jac_weighted = jac_surf * w_ip;

                    % Penalty force from contact gap difference
                    gap_current = sec_xyz_current2 - prim_xyz_current2;
                    gap_ref     = ip_map0(ip).gap;
                    gapdiff     = gap_current - gap_ref;

                    edata2.all_ip_map{1, 1}{iface, 1}(ip).gap_current = norm(gapdiff);

                    T_ip = eps * gapdiff;

                    % Debug/contact information
                    gap_info = [gap_info; iface, ip, ...
                        gap_current(1),      gap_current(2),      gap_current(3), ...
                        sec_xyz_current2(1), sec_xyz_current2(2), sec_xyz_current2(3), ...
                        prim_xyz_current2(1), prim_xyz_current2(2), prim_xyz_current2(3)];

                    % Tie-contact virtual work contribution
                    surface_TieVW = surface_TieVW + dot(T_ip, delta_g) * jac_weighted;
                    surface_Area  = surface_Area  + jac_weighted;

                end

                TieVW = TieVW + surface_TieVW;

                % Store face contribution in global vector
                TieVW_vector(iface, param_ind) = surface_TieVW;

                TieArea = TieArea + surface_Area;

            end

        end

        TieVW_out(param_ind) = TieVW;

    end


    %% --------------------------------------------------------------------
    %  Nodal force Virtual Work
    % ---------------------------------------------------------------------
    % Computes virtual work from nodal forces on interface/traction subsets.

    if hasNodalCache

        nodal_VW = nodal_VW_out(param_ind);

    else

        % Number of traction surfaces/interfaces
        nLoad_nodal = length(edata2);

        % Virtual work contribution from each nodal-force subset
        nodal_VW_vector = zeros(nLoad_nodal, 1);

        parfor nLoad_idx = 1:nLoad_nodal

            % Surface/interface connectivity, node IDs in columns 2 to 5
            conn = edata2(nLoad_idx).surface_traction.conn(:, 2:5);

            % Unique nodes on this interface
            interface_nodes = unique(conn(:));

            % Keep only valid positive node IDs
            interface_nodes = interface_nodes(interface_nodes > 0);

            % Sum nodal_force dot virtual_displacement over this node subset
            work_subset = sum(sum( ...
                nodal_force(interface_nodes, :) .* delu(interface_nodes, :)));

            nodal_VW_vector(nLoad_idx) = work_subset;

        end

        % Total nodal virtual work
        nodal_VW = sum(nodal_VW_vector);

        nodal_VW_out(param_ind) = nodal_VW;

    end

    


    %% --------------------------------------------------------------------
    %  Cost function
    % ---------------------------------------------------------------------
    % Virtual work residual:
    % internal work should balance all external contributions.
    residual = IVW - EVW - TieVW - nodal_VW;

    % Normalized squared residual contribution for this parameter set
    add_cost = (residual / (EVW+10^-8))^2;

    % Accumulate over parameter sweep
    cost_func = cost_func + add_cost;

    fprintf('EVW/IVW calculation took %.4f s\n', toc);

end


%% ------------------------------------------------------------------------
%  Final cost
% -------------------------------------------------------------------------

% Final L2-like cost over all parameter sets
cost_func = sqrt(cost_func);


%% ------------------------------------------------------------------------
%  Write FEBio simulation caches if newly computed
% -------------------------------------------------------------------------

if NewVirtualWorkFlag == 1

    tic
    writematrix(nodedat_out, nodeCachePath);
    writematrix(elemdat_out, elemCachePath);
    fprintf('writing nodedat/elemdat cache took %.4f s\n', toc);

end


%% ------------------------------------------------------------------------
%  Write EVW/TieVW/Nodal caches
% -------------------------------------------------------------------------

if ~isfile(energyFile) || ForwardCount == 1
    writematrix(EVW_out, energyFile);
end

if ~isfile(tieFile) || ForwardCount == 1
    writematrix(TieVW_out, tieFile);
end

if ~isfile(nodalFile) || ForwardCount == 1
    writematrix(nodal_VW_out, nodalFile);
end

end