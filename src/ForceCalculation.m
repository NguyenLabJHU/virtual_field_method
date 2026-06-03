classdef ForceCalculation
    methods(Static)
        function W = build_internal_force_matrix(model, edata, gauss_order, nDim)
            % elements: (nElem, 9), ecoords: (nNodes, 3), gauss_points: (nGauss, 3), weights: (nGauss, 1)
            step =  length(edata.steps);
            if isfield(model, 'elements_removed') && ...
               ~isempty(model.elements_removed) && ...
               iscell(model.elements_removed)
                elements = model.elements_removed{1};
            else
                elements = model.elements;
            end
            ecoords = edata.steps{1,step}.results.ecoords_prestress;
            [gauss_points, weights] = get_gauss_points(nDim,gauss_order);
            
            nElem = size(elements, 1);
            nGauss = size(gauss_points, 1);
            nNodes = size(ecoords, 1);
            
            reorder = [1 2 3 6 4 5]; % MATLAB indices (1-based)
            rows = [];
            cols = [];
            vals = [];
            
            for elem = 1:nElem
                node_ids = elements(elem, 2:end); % Já está em 1-based
                elem_coords = ecoords(node_ids, :); % (8, 3)
                for gp = 1:nGauss
                    xi   = gauss_points(gp, 1);
                    eta  = gauss_points(gp, 2);
                    zeta = gauss_points(gp, 3);
                    w_gp = weights(gp);

                    xi_eta_zeta = [xi eta zeta];

                    [~, dN_dxi, dN_deta, dN_dzeta] = hex8_shape_functions(xi_eta_zeta);
                
                    dX_dxi   = dN_dxi' * elem_coords;
                    dX_deta  = dN_deta' * elem_coords;
                    dX_dzeta = dN_dzeta' * elem_coords;
                    jacmat = [dX_dxi(:)'; dX_deta(:)'; dX_dzeta(:)']'; % (3,3)
                    detJ = det(jacmat);
                    [Ji,~] = compute_inv_transpose_jacobian(jacmat);
        
                    for i = 1:8
                        Gr = dN_dxi(i);
                        Gs = dN_deta(i);
                        Gt = dN_dzeta(i);
                        Gx = Ji(1,1)*Gr + Ji(1,2)*Gs + Ji(1,3)*Gt;
                        Gy = Ji(2,1)*Gr + Ji(2,2)*Gs + Ji(2,3)*Gt;
                        Gz = Ji(3,1)*Gr + Ji(3,2)*Gs + Ji(3,3)*Gt;
        
                        n = node_ids(i);
                        col_base = (elem-1)*nGauss*6 + (gp-1)*6;
        
                        % f_x: sigma_xx*Gx + sigma_xy*Gy + sigma_zx*Gz
                        rows = [rows; repmat(n*3-2,3,1)];
                        cols = [cols; col_base+ (1); col_base+(4); col_base+(6)];
                        vals = [vals; -detJ*w_gp*Gx; -detJ*w_gp*Gy; -detJ*w_gp*Gz];
        
                        % f_y: sigma_yy*Gy + sigma_xy*Gx + sigma_yz*Gz
                        rows = [rows; repmat(n*3-1,3,1)];
                        cols = [cols; col_base+(2); col_base+(4); col_base+(5)];
                        vals = [vals; -detJ*w_gp*Gy; -detJ*w_gp*Gx; -detJ*w_gp*Gz];
        
                        % f_z: sigma_zz*Gz + sigma_yz*Gy + sigma_zx*Gx
                        rows = [rows; repmat(n*3,3,1)];
                        cols = [cols; col_base+(3); col_base+(5); col_base+(6)];
                        vals = [vals; -detJ*w_gp*Gz; -detJ*w_gp*Gy; -detJ*w_gp*Gx];
                    end
                end
            end
            n_stress = nElem * nGauss * 6;
            W = sparse(rows, cols, vals, nNodes*3, n_stress);
        end



        function nodal_forces = apply_internal_force_matrix(edata, W)
            % stress_gp: (nElem, nGauss, nStress)
            % W: sparse (nNodes*3, nElem*nGauss*6)

            steps = edata.steps;
            last_step = length(steps);
            stress_gp = edata.steps{1,last_step}.results.stress_prestress;
            reorder = [1 2 3 6 4 5];
            stress_gp_perm = stress_gp(:,:,reorder);         % reorganiza componentes
            stress_gp_reshape = permute(stress_gp_perm, [3 2 1]);     % comp, gp, elem
            stress_gp_flat = reshape(stress_gp_reshape, [], 1);       % vetor coluna, igual Python flatten
            forces_flat = W * stress_gp_flat;
            nNodes = length(forces_flat)/3;
            nodal_forces = reshape(forces_flat, 3, nNodes)';
        end

        function f_boundary = compute_surface_pressure_forces_exp(model, edata, p_app, gauss_order)
            % model.surfacesp - cell array, each entry is [element_number, node1, node2, node3, node4]
            % ecoords - [nNodes x 3]
            % pressure_value - [2 x nPressureSamples] (primeira linha: pressão; segunda linha: tempo)
            % gauss_order - quadratura (ex: 2)
            % time - instante          
            
            if isfield(model, 'surfacesp_removed') && ...
               ~isempty(model.surfacesp_removed) && ...
               iscell(model.surfacesp_removed)
                surfacesp = model.surfacesp_removed{1};
            else
                surfacesp = model.surfacesp{1}; % [nFaces x 5], cada linha: [element_number, node1, node2, node3, node4]
            end
            
            % Extração dos vetores de tempo e pressão
            pressure_at_t = p_app(1);
            step =  length(edata.steps);
            ecoords = edata.steps{1,step}.results.ecoords_prestress;
        
            % Interpolação linear da pressão
            % pressure_at_t = interp1(time_vec, pressure_vec, time, 'linear', 'extrap');
        
            % Pré-alocação da matriz de forças
            nNodes = size(ecoords,1);
            f_boundary = zeros(nNodes,3);
        
            % Pontos de Gauss e pesos na face quadrilátera
            [gauss_pts, gauss_weights] = get_gauss_points(2,gauss_order); % [ngp x 2], [ngp x 1]
            nGauss = length(gauss_weights);
        
            % Carregando as faces onde a pressão é aplicada
            num_faces = size(surfacesp,1);
        
            for f = 1:num_faces
                surf_data = surfacesp(f,:);
                node_ids = surf_data(2:5); % 1-based já!
                face_coords = ecoords(node_ids, :); % [4 x 3]
 
        
                for gp = 1:nGauss
                    xi_eta = gauss_pts(gp,:); % local da face [1 x 2]
                    [N_face, dN] = evaluateShapeFunctions(xi_eta); % N_face: [4 x 1]
                    dN_dxi = dN{1};
                    dN_deta = dN{2};

                    % Jacobiano da face no ponto de Gauss
                    dx_dxi = dN_dxi' * face_coords;   % [1 x 3]
                    dx_deta = dN_deta' * face_coords; % [1 x 3]
        
                    normal_at_gp = (cross(dx_dxi, dx_deta));
            
                    % Pressão neste ponto
                    p_gp = pressure_at_t;
                    w = gauss_weights(gp);
        
                    % Acumular força nodal
                    for loc_node = 1:4
                        global_node = node_ids(loc_node); % 1-based
                        
                        % f_gp = shape * pressão * normal * peso
                        f_gp = N_face(loc_node) * p_gp * normal_at_gp * w;
        
                        f_boundary(global_node,:) = f_boundary(global_node,:) + f_gp;
                    end
                end
            end
        end

        function edata2 = overallStress (model, edata2,gauss_order_hex,nDim,matparam,prestress_time)   
            
            if isfield(model, 'elements_removed') && ...
               ~isempty(model.elements_removed) && ...
               iscell(model.elements_removed)
                elements = model.elements_removed{1};
                elemmat = model.elemmat_removed{1};
                all_elements = model.all_elements;
            else
                elements = model.elements;
                elemmat = model.elemmat;
                all_elements = model.elements;
            end

            nel  = length(elements);           % Number of volume elements
            
            % Time step array from experimental/FE results
            times           = edata2.times;
            last_time_step  = length(times);         % Use last time step (final configuration)

            ecoords = edata2.steps{1,last_time_step}.results.ecoords_prestress;      % Node coords at tn
               
            % Get gauus point data
            [gauss_points_elem, ~] = get_gauss_points(nDim, gauss_order_hex);

            % Retrieve cumulative Fpre from experimental data at final time step
            F_all = edata2.steps{1,last_time_step}.results.eFel_calc;
            %[F_all, ~] = compute_deformation_gradient(model, ecoords0, edisp, gauss_order_hex);
            
            % Retrieve cumulative Fpre from experimental data at final time step
            if prestress_time>10^-6
               F_all_pre = edata2.steps{1,last_time_step}.results.eFpre_calc; % shape: [nel x nGauss x 3 x 3]
            end
                
            % Find the index of the subset of elements wrt to the full set
            [~, map] = ismember(elements(:,1), all_elements(:,1));

            for k = 1:nel
                ielem = elements(k,1);
                global_idx = map(k);
                imatprop = zeros(size(model.matprop,2),1); 
                
                % Assemble each element's 8 node coordinates
                X = zeros(8,3);
                for node_idx = 1:8
                    X(node_idx,:) = ecoords(elements(k, node_idx+1), :);
                end
        
                % Material assignment for element
                matnum = elemmat(k);
                imatprop(:) = matparam(matnum,:);
        
        
                % ---- Loop over all Gauss points ----
                for gp = 1:size(gauss_points_elem, 1)
        
                    % Extract deformation gradient and prestrain at this element/Gauss pt
                    if abs(prestress_time) > 1e-6
                        mFbar = squeeze(F_all(global_idx,gp,:,:));           % from simulation
                        mFpre = squeeze(F_all_pre(global_idx,gp,:,:));       % from experiment
                        F     = mFbar * mFpre;                      % total F
                    else
                        mFbar = squeeze(F_all(k,gp,:,:));           % from simulation
                        mFpre = eye(3);                             % from experiment
                        F     = mFbar * mFpre;                      % total F
                    end    
        
                    % Cauchy stress tensor
                    sigma =  computeStress(F, model,matnum, imatprop,ielem,gp);
                    edata2.steps{1,last_time_step}.results.stress_prestress(k,gp,:) = [sigma(1,1); sigma(2,2); sigma(3,3); ...
                                sigma(2,3); sigma(1,3); sigma(1,2)];
                end
            end
        end
    end
end


function [Ji, detJ] = compute_inv_transpose_jacobian(jacmat)
% COMPUTE_INV_TRANSPOSE_JACOBIAN Calculates the inverse transpose of a 3x3 Jacobian matrix analytically.
% 
% This approach avoids the computational overhead and floating-point errors
% of the built-in inv() function, which is crucial for stability and speed 
% in Finite Element Analysis (FEA) integration loops.
%
% Inputs:
%   jacmat - 3x3 Jacobian matrix evaluated at a Gauss integration point
%
% Outputs:
%   Ji   - 3x3 Inverse transpose of the Jacobian matrix (J^-T)
%   detJ - Determinant of the Jacobian matrix (used for volume integration)

    % 1. Extract matrix components directly to avoid repeated indexing overhead
    J11 = jacmat(1,1); J12 = jacmat(1,2); J13 = jacmat(1,3);
    J21 = jacmat(2,1); J22 = jacmat(2,2); J23 = jacmat(2,3);
    J31 = jacmat(3,1); J32 = jacmat(3,2); J33 = jacmat(3,3);

    % 2. Analytical calculation of the determinant
    detJ = J11*(J22*J33 - J23*J32) - J12*(J21*J33 - J23*J31) + J13*(J21*J32 - J22*J31);

    % 3. Check for negative or near-zero determinant (element inversion or collapse)
    % This is a crucial safeguard to prevent silent errors during Newton-Raphson iterations
    if detJ <= 1e-12
        error('Negative or zero Jacobian determinant (detJ = %g). Element is highly distorted or inverted.', detJ);
    end

    % 4. Calculate Cofactor Matrix components (C)
    C11 =  (J22*J33 - J23*J32);
    C12 = -(J21*J33 - J23*J31);
    C13 =  (J21*J32 - J22*J31);

    C21 = -(J12*J33 - J13*J32);
    C22 =  (J11*J33 - J13*J31);
    C23 = -(J11*J32 - J12*J31);

    C31 =  (J12*J23 - J13*J22);
    C32 = -(J11*J23 - J13*J21);
    C33 =  (J11*J22 - J12*J21);

    % 5. Final assembly of the Inverse Transpose (Ji = J^-T)
    % Note: The inverse transpose is exactly the cofactor matrix divided by the determinant.
    % We do not need to transpose C here because the math naturally aligns.
    Ji = (1 / detJ) * [C11, C12, C13; 
                       C21, C22, C23; 
                       C31, C32, C33];
end







