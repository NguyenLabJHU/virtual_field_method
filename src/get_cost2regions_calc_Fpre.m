function cost = get_cost2regions_calc_Fpre(path, mymodel, model, simp_model,edata,... 
    x, p_app,gauss_order,prestress_time,eps,changing_matrix,Normalizer, ...
    ops_matrix_struct, run_simple_model, mat_surface_traction,rParts,nDim)
    
    global totalRunCount edata_with_Fpre_step ForwardCount edata_with_Fpre_step_simp
     


    %Creating material parameters vector to sweep on the PVW calculation
    [ground_truth_mat,matparam_sweep,matparam_complete] = SweepMatrix(model,...
        changing_matrix,x,Normalizer,ops_matrix_struct);


    % % Kick off parallel pool
    % if isempty(gcp('nocreate'))
    % 
    %     N = str2double(getenv('SLURM_CPUS_PER_TASK'));
    % 
    %     if isnan(N) || N == 0
    %         N = feature('numcores'); % fallback se rodar local
    %     end
    % 
    %     parpool('local', N);
    % end

    
    % Try-catch block: handle potential failures in FEBio or cost calculation routines.
    try
        % Accumulate prestress calculation for given material parameters updated at every set of parameter.
        % Only update Fpre every 10 evaluations of the cost function
        if ForwardCount == 1
            % Always run on the first call, and then every time a minimun
            % is reached
            mydir_data = path.data;
            
            % Calculating the prestress
            edata_pre = accumulate_Fpre_from_edata(mydir_data, mymodel, ...
                gauss_order, prestress_time, matparam_complete, edata,model);
            
            % Save the new edata
            edata_with_Fpre_step = edata_pre;

            % Calculate the simplified version traction
            if run_simple_model == "True"
                    
              % --- Find the matrix to extract the nodal force from each element ---
              W_matrix = ForceCalculation.build_internal_force_matrix(simp_model,edata_pre, gauss_order,nDim);
              edata_simp = ForceCalculation.overallStress(simp_model, edata_pre,gauss_order,nDim,matparam_complete,prestress_time);
              nodal_forces = ForceCalculation.apply_internal_force_matrix(edata_simp, W_matrix);
              f_boundary = ForceCalculation.compute_surface_pressure_forces_exp(simp_model, edata_simp, p_app, gauss_order);       
              nodal_forces = nodal_forces + f_boundary;

               % W_matrix = ForceCalculation.build_internal_force_matrix(model,edata_pre, gauss_order,nDim);
               % edata_pre = ForceCalculation.overallStress(model, edata_pre,gauss_order,nDim,matparam_complete,prestress_time);
               % nodal_forces = ForceCalculation.apply_internal_force_matrix(edata_pre, W_matrix);
               % f_boundary = ForceCalculation.compute_surface_pressure_forces_exp(model, edata_pre, p_app, gauss_order);
               %      % 
               % nodal_foces_net = nodal_forces + f_boundary;

                edata_simp.nodal_forces = nodal_forces;
                edata_with_Fpre_step_simp = edata_simp;
            end

            %edata_change = updateEdataPrestress(model, edata, edata_pre, rParts);
            
            edata = edata_pre;

            ForwardCount = 1;
            
            % Uncomment for debug:
            fprintf('[Fpre updated at call %d]\n', totalRunCount);
            fprintf('[Number of Forward count call %d]\n', ForwardCount);
            
        else
            % Use the previous Fpre to save on calculation
            edata = edata_with_Fpre_step;

            edata_simp = edata_with_Fpre_step_simp;

            nodal_forces = edata_simp.nodal_forces;
            
            % Uncomment for debug:
            fprintf('[Fpre reused at call %d]\n', totalRunCount);
        end


       if run_simple_model == "True" 
            % Compute the cost function using virtual work integrals for the
            % simplified model      
            [~, ~, ~, cost_func] = calc_virtual_work_variation_integration_simp(path, mymodel, simp_model, ...
                edata_simp, matparam_complete, matparam_sweep ,ground_truth_mat,p_app,...
                gauss_order,eps,changing_matrix,ops_matrix_struct,model,nodal_forces);

       else
            % Compute the cost function using virtual work integrals.
            [~, ~, ~, cost_func] = calc_virtual_work_variation_integration2(path, mymodel, model, ...
                edata, matparam_complete, matparam_sweep ,ground_truth_mat,p_app,...
                gauss_order,eps,changing_matrix,ops_matrix_struct);
       end

        
        % Check result validity: cost must be finite and non-NaN.
        if isnan(cost_func) || ~isfinite(cost_func)
            error('NaN or Inf encountered in cost function. Triggering restart.');
        end

        % Take logarithm of cost (for scale, stability, or easier optimization).
        cost = log10(cost_func)

        % Debug print (only when totalRunCount == 1 for initialization/testing).
        if totalRunCount == 1
            fprintf('total count %.4f s\n', totalRunCount);
        end

        % Increment global run counter (tracks number of cost calls/runs).
        totalRunCount = totalRunCount+1;
        ForwardCount = ForwardCount+1;

    catch ME
        % Print error message for debugging if something fails in try block.
        fprintf('Restarting due to error: %s\n', ME.message);
        % Rethrow the error so MultiStart/fmincon can handle and possibly restart.
        rethrow(ME);  
    end




    
end