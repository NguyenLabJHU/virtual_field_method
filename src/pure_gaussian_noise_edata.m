function edata_noisy = pure_gaussian_noise_edata(edata, model, noise_percent, sigma_additive)
    % pure_gaussian_noise_checked: Applies pure white Gaussian noise
    % and checks the Jacobian to ensure no elements invert.
    
    edata_noisy = edata;
    nsteps = size(edata.steps, 2);
    coords_ref = edata.steps{1, 1}.results.ecoords;
    [num_nodes, dim] = size(coords_ref);
    
    % Extract connectivity (assuming hex)
    conn = model.elements(:, 2:9);

    % --- 1. Defining Displacement-Based Scale ---
    all_disps = [];
    for s = 1:nsteps
        all_disps = [all_disps; edata.steps{s}.results.edisp]; 
    end
    max_u = max(abs(all_disps(:)));
    noise_std = (noise_percent/100) * max_u + sigma_additive;

    % --- 2. Noise Injection ---
    fprintf('Starting pure Gaussian noise injection with Jacobian check...\n');
    for step_idx = 1:nsteps
        ecoords_orig = edata.steps{1, step_idx}.results.ecoords;
        
        success = false;
        attempts = 0;
        current_scale = 1.0; % Starts with 100% of the requested noise scale

        % Attempts up to 10 times to reduce noise if elements invert
        while ~success && attempts < 10
            % Generates PURE white noise (independent for each node)
            noise_raw = (noise_std * current_scale) .* randn(num_nodes, dim);
            
            % Applies to coordinates 
            test_coords = ecoords_orig + noise_raw;
            
            % --- JACOBIAN CHECK ---
            if check_jacobian(test_coords, conn)
                % Passed the test! Accepts coordinates and exits the loop
                edata_noisy.steps{1, step_idx}.results.ecoords = test_coords;
                success = true;
            else
                % Failed.
                attempts = attempts + 1;
                current_scale = current_scale * 0.5;
                fprintf('Warning (Step %d): Element inverted. Reducing noise scale to %.4f\n', step_idx, current_scale);
            end
        end
        
        % If the mesh breaks even after 10 attempts, keeps the original data
        if ~success
            warning('Step %d: Impossible to inject noise without inverting the mesh. Keeping original coordinates.', step_idx);
            edata_noisy.steps{1, step_idx}.results.ecoords = ecoords_orig;
        end
    end

    % --- 3. Final Kinematic Consistency ---
    % Updates displacements (u = x - X) to match the noisy coordinates
    coords1_noisy = edata_noisy.steps{1, 1}.results.ecoords;
    for step_idx = 1:nsteps
        curr_noisy = edata_noisy.steps{1, step_idx}.results.ecoords;
        edata_noisy.steps{1, step_idx}.results.edisp = curr_noisy - coords1_noisy;
    end
    
    fprintf('Noise injection completed.\n');
end

% =========================================================================
% LOCAL FUNCTION: JACOBIAN CHECK
% =========================================================================
function elements_ok = check_jacobian(coords, conn)
    % coords: matrix [nodes x 3] with current coordinates (x, y, z)
    % conn:   matrix [elements x 8] with mesh connectivities
    
    num_elements = size(conn, 1);
    elements_ok = true; 
    
    % Shape function derivatives at the center of the hexahedral element (C3D8)
    dN = 0.125 * [
        -1, -1, -1;
         1, -1, -1;
         1,  1, -1;
        -1,  1, -1;
        -1, -1,  1;
         1, -1,  1;
         1,  1,  1;
        -1,  1,  1
    ];

    for i = 1:num_elements
        node_indices = conn(i, :);
        element_coords = coords(node_indices, :); 
        
        % Jacobian Matrix at the center of the element
        J = dN' * element_coords;
        
        % Checks the determinant
        if det(J) <= 0
            elements_ok = false;
            return; % Exits on the first error found to optimize speed
        end
    end
end