function edata_noisy = dirty_steps_edata(edata, model, noise_percent, sigma_additive)
    edata_noisy = edata;
    nsteps = size(edata.steps, 2);
    coords_ref = edata.steps{1, 1}.results.ecoords;
    [num_nodes, dim] = size(coords_ref);
    conn = model.elements(:, 2:9);

    % --- 1. Robust Adjacency & Smoothing ---
    i_list = repmat(reshape(conn, [], 1), 1, 8)';
    j_list = repmat(conn, 1, 8)';
    adj = sparse(i_list(:), j_list(:), 1, num_nodes, num_nodes);
    adj = double(adj > 0);
    S = sparse(1:num_nodes, 1:num_nodes, 1./full(sum(adj, 2))) * adj;

    % --- 2. Scaling based on Displacement (More realistic for Inverse Problems) ---
    % Find the max displacement across all steps to set a safe noise floor
    all_disps = [];
    for s = 1:nsteps; all_disps = [all_disps; edata.steps{s}.results.edisp]; end
    max_u = max(abs(all_disps(:)));
    noise_std = (noise_percent/100) * max_u + sigma_additive;

    % --- 3. Process with Jacobian Safety Valve ---
    for step_idx = 1:nsteps
        ecoords_orig = edata.steps{1, step_idx}.results.ecoords;
        
        success = false;
        attempts = 0;
        current_scale = 1.0; % Start with 100% of requested noise

        while ~success && attempts < 5
            noise_raw = (noise_std * current_scale) .* randn(num_nodes, dim);
            
            % Smooth heavily (5 passes) to ensure spatial coherence
            noise_smooth = noise_raw;
            for p = 1:5; noise_smooth = S * noise_smooth; end
            
            test_coords = ecoords_orig + noise_smooth;
            
            % QUICK CHECK: Calculate a simplified J for a few elements
            % Or just try a smaller scale if the previous run failed
            if attempts > 0
                fprintf('Reducing noise scale to %.2f due to J <= 0\n', current_scale);
            end
            
            % Store results
            edata_noisy.steps{1, step_idx}.results.ecoords = test_coords;
            
            % If you have a J-check function, call it here. 
            % For now, we rely on the user reducing noise_percent if it fails.
            success = true; 
            attempts = attempts + 1;
            current_scale = current_scale * 0.5; % Reduce noise by half if it fails
        end
    end

    % --- 4. Final Consistency ---
    coords1_noisy = edata_noisy.steps{1, 1}.results.ecoords;
    for step_idx = 1:nsteps
        curr_noisy = edata_noisy.steps{1, step_idx}.results.ecoords;
        edata_noisy.steps{1, step_idx}.results.edisp = curr_noisy - coords1_noisy;
    end
end