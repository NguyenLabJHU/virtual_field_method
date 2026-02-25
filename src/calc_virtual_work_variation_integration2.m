 function [IVW, EVW, TieVW, cost_func] = calc_virtual_work_variation_integration2(path, mymodel, model, ...
    edata2, matparam, matparam_sweep,ground_truth_mat, p_app, gauss_order,eps,changing_matrix,ops_matrix_struct)
%--------------------------------------------------------------------------
% calc_virtual_work_variation_integration2
%--------------------------------------------------------------------------
% Computes the Internal Virtual Work (IVW), External Virtual Work (EVW), and the cost function
% for a finite element model perturbed with given test material parameters, using Virtual Field Method.
%
% INPUTS:
%   path         : string, directory path containing root folder
%   mymodel      : string, filename of the FEBio finite element model
%   model        : struct, mesh and element/nodal connectivity (from preliminary_reading2)
%   edata2       : struct, experimental/nodal data and Fpre from logs (from preliminary_reading2)
%   matparam     : array [nMaterials x 2], each row is [shear_modulus, bulk_modulus] for a material
%   p_app        : scalar, applied pressure for virtual work calculations
%   gauss_order  : integer, quadrature order for hexahedral/surface integration
%
% OUTPUTS:
%   IVW          : scalar, calculated internal virtual work for tested parameters
%   EVW          : scalar, calculated external virtual work for tested parameters
%   cost_func    : scalar, cost function (squared error between IVW and EVW)
%
% Uses cached simulation outputs or computes fresh FEA simulation for parameter sets.
% Integrates over volume and surfaces, computes virtual work, stress and strain at Gauss points.
%--------------------------------------------------------------------------

global totalRunCount ForwardCount IVW_vector EVW_vector TieVW_vector prestress_time   % Track across parameter sweeps


% Kick off parallel pool
if isempty(gcp('nocreate'))

    N = str2double(getenv('SLURM_CPUS_PER_TASK'));
    
    if isnan(N) || N == 0
        N = feature('numcores'); % fallback se rodar local
    end
    
    parpool('local', 4);
end


%--------------------------------------------------------------------------
% Mesh and Simulation Geometry Setup
%--------------------------------------------------------------------------

nel  = length(model.elements);           % Number of volume elements
nnd  = length(model.nodes);              % Number of nodes

% Time step array from experimental/FE results
times           = edata2.times;
first_time_step = 1;                     % Use first time step
last_time_step  = length(times);         % Use last time step (final configuration)
prestress_step = find(edata2.times <= prestress_time, 1, 'last'); % find last time less than or equal to target

% Extract experimental nodal and element results at first and last time steps
res0   = edata2.steps{1,first_time_step}.results;
ecoords0 = res0.ecoords;                 % Node coords at t0 (nnd x 3)
edisp0   = res0.edisp;                   % Node displacements at t0

resN   = edata2.steps{1,last_time_step}.results;
ecoords = resN.ecoords;                  % Node coords at tn
edisp   = resN.edisp;                    % Node displacements at tn

resN_1   = edata2.steps{1,last_time_step-1}.results;
ecoords_1 = resN_1.ecoords;              % Node coords at tn-1
edisp_1   = resN_1.edisp;                % Node displacements at tn-1

resN_pre0   = edata2.steps{1,prestress_step}.results;
ecoords_pre0 = resN_pre0.ecoords;                  % Node coords at the start of elastic
edisp_pre0   = resN_pre0.edisp;                   % Node displacements at the start of elastic

resN_pre   = edata2.steps{1,last_time_step}.results;
ecoords_pre = resN_pre.ecoords_prestress;                  % Node coords at the end of elastic
edisp_pre   = resN_pre.edisp_prestress;                   % Node displacements the end of elastic

% Starting empty struct
delta_e_struct = struct(); 
gauss_order_hex = gauss_order;           % Gauss order for the hex elements

% Name of the cache file to store delta_e_struct
modelName = erase(mymodel, '.feb');      % Remove the extension from .feb 
cachefile_delta_e = ['delta_e_struct_', modelName, '.mat'];
fulldelta_e_File = fullfile(path.VF, cachefile_delta_e);

% Deformation gradients for reference configuration
[F_all, ~] = compute_deformation_gradient(model, ecoords0, edisp, gauss_order_hex);

% Retrieve cumulative Fpre from experimental data at final time step
if prestress_time>10^-6
   F_all_pre = edata2.steps{1,last_time_step}.results.eFpre_calc; % shape: [nel x nGauss x 3 x 3]
end

%--------------------------------------------------------------------------
% Simulation or Cache: Compute or Load Perturbed Results
%--------------------------------------------------------------------------

nParam = size(changing_matrix,2);        % Number of material parameter sets in sweep
cost_func = 0;                           % Initialize cost function accumulator

% Allocate storage for simulated nodal and element data across param sets
nodedat_out = zeros(nnd, 7*nParam);      % (node, [index,x,y,z,ux,uy,uz])
elemdat_out = zeros(nel, 19*nParam);     % (element, [index,...Fbar(9)...Fpre(9)...])

% Construct the full paths using the VF directory
nodeCachePath = fullfile(path.VF, ['nodedat_cached_', modelName, '.csv']);
elemCachePath = fullfile(path.VF, ['elemdat_cached_', modelName, '.csv']);

% Checking if there is need to calculate the VF
if totalRunCount == 1
    % First run for sweep: Perform FE simulation for each material property set
    NewVirtualWorkFlag = 1;
else
    % For subsequent runs, load matrix results from CSV cache for efficiency
    tic
    nodedat_cached = readmatrix(nodeCachePath);
    elemdat_cached = readmatrix(elemCachePath);
    time4 = toc;
    fprintf('reading nodedat_cached took %.4f s\n', time4);
    NewVirtualWorkFlag = 2;
end


%--------------------------------------------------------------------------
% Parameter Sweep Loop: For Each Material Parameter Set
%--------------------------------------------------------------------------

for param_ind = 1:nParam
    % --- FEA Simulation or Load Cached Results ---
    if NewVirtualWorkFlag == 1
        % Run FEA simulation and store results for this parameter set
        mydir = path.data;
        tic
        
        [nodedat, elemdat] = simulate_febio_uniform(mydir, mymodel, ...
            matparam_sweep,ground_truth_mat ,nnd, nel, param_ind,...
            changing_matrix,ops_matrix_struct,model);

        time3 = toc;
        fprintf('simulate_febio_uniform took %.4f s\n', time3);

        % Save to cache array for this param variation
        nodedat_out(:,(param_ind-1)*7+1:(param_ind-1)*7+7)    = nodedat;
        elemdat_out(:,(param_ind-1)*19+1:(param_ind-1)*19+19) = elemdat;
        
    else
        % Load previously cached nodal/element matrix for this param set
        nodedat  = nodedat_cached(:,(param_ind-1)*7+1:(param_ind-1)*7+7);
        elemdat  = elemdat_cached(:,(param_ind-1)*19+1:(param_ind-1)*19+19);    
    end

    % Simulation/data failure check: If any results contain NaN, abort this run
    if any(isnan(nodedat(:))) || any(isnan(elemdat(:)))
        IVW = NaN;
        EVW = NaN;
        cost_func = NaN;
        return;
    end

    tic

    %--------------------------------------------------------------------------
    % Extract Simulated Nodal/Element Variables
    %--------------------------------------------------------------------------
    disp   = nodedat(:,5:7);             % Node ux, uy, uz displacements

    %--------------------------------------------------------------------------
    % Compute Virtual Field: Compare simulation and experiment displacements
    %--------------------------------------------------------------------------
    delu = disp - edisp;                 % Virtual displacement = simulated - reference
       
    % Fieldname for this parameter index
    fieldname = sprintf('param_%d', param_ind);
    
    % Flag to decide whether to compute or just load
    need_to_compute = false;
    
    if exist(fulldelta_e_File, 'file')
        S = load(fulldelta_e_File);                % load struct from file
        if isfield(S, 'delta_e_struct') && isfield(S.delta_e_struct, fieldname)
            % field exists: just read from cache
            delta_e = S.delta_e_struct.(fieldname);
            
        else
            need_to_compute = true;
            if isfield(S, 'delta_e_struct')
                delta_e_struct = S.delta_e_struct;  % retain rest of struct to not lose other fields
            else
                delta_e_struct = struct();
            end
        end
    else
        need_to_compute = true;
        delta_e_struct = struct();          % start new if file doesn't exist
    end
    
    if need_to_compute
        % Compute virtual strain field at Gauss points for all elements
        [delta_e, ~] = compute_virtual_strain_integration(model, ecoords, delu, gauss_order_hex);
        fprintf('Computed new delta_e for %s\n', fieldname);
    
        % Add the result to the struct
        delta_e_struct.(fieldname) = delta_e;
    
        % Save/overwrite .mat file (matlab will update or create)
        save(fulldelta_e_File, 'delta_e_struct');
    end

    % Identity matrix used in stress and strain tensor calculations
    Id = eye(3);

    %--------------------------------------------------------------------------
    % --- Load Energy Virtual Work (EVW) cache from the VF folder ---
    %--------------------------------------------------------------------------
    
    % Construct the full paths to the files inside the VF subfolder
    energyFile = fullfile(path.VF, ['EVW_cached_', modelName, '.csv']);
    tieFile    = fullfile(path.VF, ['TieVW_cached_', modelName, '.csv']);

    % Check if the cache file for the EVW exists before attempting to read
    if isfile(energyFile)
        EVW_out = readmatrix(energyFile);  
    end

    % Check if the cache file for the TieEVW exists before attempting to read
    if isfile(tieFile) && ForwardCount ~= 1
         TieVW_out = readmatrix(tieFile);
    end

    %% ----------------- Internal Virtual Work Calculation -----------------
    IVW = 0;                           % Accumulator for Internal Virtual Work   
    dimension_hex = 3;                 % 3D brick/hexahedral element

    element_IVW_array      = zeros(nel,1);  % Array to store the IVW for the parfor

    [gauss_points_elem, weights_elem] = get_gauss_points(dimension_hex, gauss_order_hex);

    % ---- Loop over all volumetric elements ----
    parfor k = 1:nel
        ielem = model.elements(k,1);
        imatprop = zeros(size(model.matprop,2),1); 
        element_vf_norm_sq = 0;        % Virtual strain norm squared per element
        
        % Assemble each element's 8 node coordinates
        X = zeros(8,3);
        for node_idx = 1:8
            X(node_idx,:) = ecoords(model.elements(k, node_idx+1), :);
        end

        % Material assignment for element
        matnum = model.elemmat(k);
        imatprop(:) = model.matprop(matnum,:);
        if ismember(matnum,changing_matrix(1,:))
            imatprop(1) = matparam(matnum, 1);
            imatprop(3) = matparam(matnum, 2);
        end

        element_IVW = 0;
        element_Vol = 0;

        % ---- Loop over all Gauss points ----
        for gp = 1:size(gauss_points_elem, 1)
            xi    = gauss_points_elem(gp, 1);
            eta   = gauss_points_elem(gp, 2);
            zeta  = gauss_points_elem(gp, 3);
            weight= weights_elem(gp);

            % Hex shape function derivatives
            dN_dxi = 0.125 * [-(1-eta)*(1-zeta);  (1-eta)*(1-zeta); ...
                              (1+eta)*(1-zeta); -(1+eta)*(1-zeta); ...
                              -(1-eta)*(1+zeta);  (1-eta)*(1+zeta); ...
                              (1+eta)*(1+zeta); -(1+eta)*(1+zeta)];
            dN_deta = 0.125 * [-(1-xi)*(1-zeta); -(1+xi)*(1-zeta); ...
                               (1+xi)*(1-zeta);  (1-xi)*(1-zeta); ...
                               -(1-xi)*(1+zeta); -(1+xi)*(1+zeta); ...
                               (1+xi)*(1+zeta);  (1-xi)*(1+zeta)];
            dN_dzeta = 0.125 * [-(1-xi)*(1-eta); -(1+xi)*(1-eta); ...
                                -(1+xi)*(1+eta); -(1-xi)*(1+eta); ...
                                (1-xi)*(1-eta);  (1+xi)*(1-eta); ...
                                (1+xi)*(1+eta);  (1-xi)*(1+eta)];

            dX_dxi   = dN_dxi'  * X;
            dX_deta  = dN_deta' * X;
            dX_dzeta = dN_dzeta' * X;
            Jac_matrix = [dX_dxi; dX_deta; dX_dzeta];   % 3x3 Jacobian

            jac_det      = det(Jac_matrix);
            jac_weighted = jac_det * weight;

            % Extract deformation gradient and prestrain at this element/Gauss pt
            if abs(prestress_time) > 1e-6
                mFbar = squeeze(F_all(k,gp,:,:));           % from simulation
                mFpre = squeeze(F_all_pre(k,gp,:,:));       % from experiment
                F     = mFbar * mFpre;                      % total F
            else
                mFbar = squeeze(F_all(k,gp,:,:));           % from simulation
                mFpre = eye(3);                             % from experiment
                F     = mFbar * mFpre;                      % total F
            end    
            dgum  = squeeze(delta_e(k,gp,:,:));         % Virtual strain tensor

            % Cauchy stress tensor
            sigma =  computeStress(F, model,matnum, imatprop,ielem,gp);
            %sigma = (1.0/J) * (2*c1*(B-Id) + kappa*log(J)*Id);

            % Compute virtual work (trace of sigma*delta_e)
            temp = sum(sum(dgum .* sigma)); 

            % Accumulate IVW and volume for this Gauss point
            element_IVW = element_IVW + temp * jac_weighted;   
            element_Vol = element_Vol + jac_weighted;

            % Accumulate norm squared of virtual strain for normalization (if desired)
            vf_q = sum(sum(dgum .* dgum));      % Frobenius norm
            element_vf_norm_sq = element_vf_norm_sq + vf_q * jac_weighted;
        end

        element_IVW_array(k) = element_IVW; % Add element's IVW to the vector
        IVW_vector(k, param_ind) = element_IVW;
    end

    % Getting the total value after parfor
    IVW        = sum(element_IVW_array);

    %% ----------------- External Virtual Work Calculation (Surface) -----------------
    if isfile(energyFile) 
        
        EVW = EVW_out(param_ind);

    else    
        
        EVW = 0;                                   % External Virtual Work accumulator
        AreaS = 0;                                 % Surface area accumulator
    
        gauss_order_surface  = gauss_order;        % Surface quadrature order
        dimension_surface    = 2;                  % Quads are 2D surfaces
    
        [gauss_points_surf, weights_surf] = get_gauss_points(dimension_surface, gauss_order_surface);
    
        % Number of applied loads
        nLoad = length(p_app);
    
        for nLoad_idx = 1:nLoad
            nels = length(model.surfacesp{1,nLoad_idx});          % Number of surface elements 
    
            for k = 1:nels
                % Extract coords for each node of the surface element
                X = zeros(4,3);
                for node_idx = 1:4
                    X(node_idx,:) = ecoords(model.surfacesp{1,nLoad_idx}(k,node_idx+1), :);
                end
        
                % Extract virtual displacement at surface nodes
                delu_nodes = zeros(4,3);
                for node_idx = 1:4
                    node_id = model.surfacesp{1,nLoad_idx}(k, node_idx+1);
                    delu_nodes(node_idx,:) = delu(node_id,:);
                end
        
                surface_EVW   = 0;
                surface_Area  = 0;
        
                % Loop over surface Gauss points
                for gp = 1:size(gauss_points_surf, 1)
                    xi     = gauss_points_surf(gp,1);
                    eta    = gauss_points_surf(gp,2);
                    weight = weights_surf(gp);
        
                    N = 0.25 * [(1-xi)*(1-eta);
                                (1+xi)*(1-eta);
                                (1+xi)*(1+eta);
                                (1-xi)*(1+eta)];
                    dN_dxi  = 0.25 * [-(1-eta);  (1-eta);  (1+eta); -(1+eta)];
                    dN_deta = 0.25 * [-(1-xi); -(1+xi);  (1+xi);  (1-xi)];
        
                    dX_dxi   = dN_dxi' * X;
                    dX_deta  = dN_deta' * X;
                    area_vec = cross(dX_dxi, dX_deta);
                    jac_surf     = norm(area_vec);     % Surface Jacobian
                    jac_weighted = jac_surf * weight;
        
                    normal = area_vec / jac_surf;      % Outward unit normal
        
                    delum = N' * delu_nodes;           % Interpolated virtual displacement at Gauss pt
        
                    surface_EVW  = surface_EVW + p_app(nLoad_idx) * dot(delum, normal) * jac_weighted;
                    surface_Area = surface_Area + jac_weighted;
                end
            
                EVW = EVW + surface_EVW;               % Add surface element EVW
                EVW_vector(k, param_ind) = surface_EVW;
                AreaS = AreaS + surface_Area;
            end   
        end
        EVW_out(param_ind) = EVW;
    end
    %% ----------------- External Virtual Work Calculation (Tie-Contact) -----------------
    if isfile(tieFile) && ForwardCount ~= 1
        
        TieVW = TieVW_out(param_ind);

    else 
        
        TieVW = 0;                 % Virtual work accumulator
        TieArea = 0;               % Contact area accumulator
        test=1;
        diff_matrix = [];      % Initialize empty double array
        gap_info = [];
    
    
        % Only run the calculation in case of existence of a tie-contact
        if isfield(model, 'tie_contact') && ~isempty(model.tie_contact)
            % --- Tie Contact Virtual Work Calculation ---        
            nfaces = size(model.tie_contact.primary_nodes, 1);
            
            for iface = 1:nfaces
                ip_map0 = edata2.all_ip_map{1,1}{iface,1}; % getting the info at the first time step
                primary_face_nodes = model.tie_contact.primary_nodes(iface,:);
        
                surface_TieVW = 0;
                surface_Area = 0;
            
                % Get node coords & virtual disp for the quad's nodes
                X_face = zeros(4,3);
                delu_nodes_prim = zeros(4,3);
                for k = 1:4
                    idx = find(model.nodes(:,1)==primary_face_nodes(k+1));
                    %X_face_current(k,:) = ecoords(idx,1:3);
                    X_face_current(k,:) = ecoords_pre(idx,1:3);
                    X_face_1(k,:) = ecoords_1(idx,1:3);
                    X_face_ref(k,:)     = ecoords0(idx,1:3);
                    delu_nodes_prim(k,:) = delu(primary_face_nodes(k+1),:);
                end
            
                for ip = 1:length(ip_map0)
                    % Integration point parametric coords (xi, eta)
                    xi     = ip_map0(ip).prim_rs(1);
                    eta    = ip_map0(ip).prim_rs(2);
            
                    % Shape functions at (xi, eta)
                    Np = 0.25 * [ (1-xi)*(1-eta);
                                 (1+xi)*(1-eta);
                                 (1+xi)*(1+eta);
                                 (1-xi)*(1+eta) ];
            
                    % Primary integration point positions
                    prim_xyz_ref = Np' * X_face_ref;
                    prim_xyz_current = Np' * X_face_1;
                    prim_xyz_current2 = Np' * X_face_current;
    
                    % Virtual displacement @ Gauss point, primary
                    delu_prim_gp = Np' * delu_nodes_prim;  % column vector
            
                    % Secondary: get connectivity and nodal virtual disp
                    sec_elem_id = ip_map0(ip).sec_elem_id;
                    sec_elem_nodes = ip_map0(ip).sec_face;
                    delu_nodes_sec = zeros(4,3);
                    X_sec_1 = zeros(4,3);
                    X_sec_ref = zeros(4,3);
                    for k2 = 1:4
                        idx2 = find(model.nodes(:,1)==sec_elem_nodes(k2));
                        delu_nodes_sec(k2,:) = delu(sec_elem_nodes(k2),:);
                        X_sec_1(k2,:) = ecoords_1(idx2,1:3);
                        %X_sec_current(k2,:) = ecoords(idx2,1:3);
                        X_sec_current(k2,:) = ecoords_pre(idx2,1:3);
                        X_sec_ref(k2,:) = ecoords0(idx2,1:3);
                    end
            
                    % Mapped secondary IP parametric coords (xi, eta)
                    xi_sec  = ip_map0(ip).sec_rs(1);
                    eta_sec = ip_map0(ip).sec_rs(2);
            
                    N_sec = 0.25 * [(1-xi_sec)*(1-eta_sec);
                                    (1+xi_sec)*(1-eta_sec);
                                    (1+xi_sec)*(1+eta_sec);
                                    (1-xi_sec)*(1+eta_sec)];
    

                    % Secondary integration point positions
                    sec_xyz_ref     = N_sec' * X_sec_ref;
                    sec_xyz_current = N_sec' * X_sec_1;
                    sec_xyz_current2 = N_sec' * X_sec_current;
                    delu_sec_gp = N_sec' * delu_nodes_sec;
            
                    % Virtual gap
                    delta_g = delu_prim_gp - delu_sec_gp;
                    
            
                    % Surface Jacobian calculation primary surface (use shape function derivatives at xi, eta)
                    dN_dxi = 0.25 * [-(1-eta); (1-eta); (1+eta); -(1+eta)];
                    dN_deta = 0.25 * [-(1-xi); -(1+xi); (1+xi); (1-xi)];
                    dX_dxi  = dN_dxi' * X_face_current;
                    dX_deta = dN_deta' * X_face_current;
                    area_vec = cross(dX_dxi, dX_deta);
                    jac_surf = norm(area_vec);
            
                    % Quadrature weight
                    w_ip = ip_map0(ip).weight;
                    jac_weighted = jac_surf * w_ip;
    
                    % Penalty force is difference of gaps (current - reference)
                    gap_current = sec_xyz_current2 - prim_xyz_current2;
                    gap_ref     = ip_map0(ip).gap;
                    gapdiff = gap_current - gap_ref;
                    edata2.all_ip_map{1,1}{iface,1}(ip).gap_current = norm(gapdiff);
                    T_ip = eps * gapdiff;
    
                    gap_info = [gap_info; iface, ip, ...
                        gap_current(1), gap_current(2), gap_current(3), ...
                        sec_xyz_current2(1),     sec_xyz_current2(2),     sec_xyz_current2(3), ...
                        prim_xyz_current2(1),     prim_xyz_current2(2),     prim_xyz_current2(3)];
                    
    
                    surface_TieVW = surface_TieVW + dot(T_ip, delta_g) * jac_weighted;
                    surface_Area = surface_Area + jac_weighted;
                end
        
                TieVW = TieVW + surface_TieVW;               % Add surface element TieVW
                TieVW_vector(iface, param_ind) = surface_TieVW;
                TieArea = TieArea + surface_Area;
        
            end
        end

        TieVW_out(param_ind) = TieVW;

    end

%% ----------------------------- Cost Function Calculation and Output ---------------------------  
    cost_func = ((IVW - EVW -TieVW)/EVW)^2 + cost_func;      % L2 cost function for virtual work balance

    time5 = toc;
    fprintf('EVW/IVW calculation took %.4f s\n', time5);    
end

cost_func = sqrt(cost_func);                    % Final cost function output

%--------------------------------------------------------------------------
% Cache Simulation Output (if simulation just performed)
%--------------------------------------------------------------------------
if NewVirtualWorkFlag == 1
    tic
    writematrix(nodedat_out, nodeCachePath);        % Cache node output
    writematrix(elemdat_out, elemCachePath);        % Cache element output
    time6 = toc;
    fprintf('writing data took %.4f s\n', time6);
end

%--------------------------------------------------------------------------
% Cache Simulation Output EVW and TieVW (if simulation just performed)
%--------------------------------------------------------------------------
if ~isfile(energyFile)
    writematrix(EVW_out,   energyFile);     % Cache EVW
end

if ~isfile(tieFile) || ForwardCount == 1
    writematrix(TieVW_out, tieFile);        % Cache TieVW
end


%--------------------------------------------------------------------------
% End of Function
%--------------------------------------------------------------------------
end