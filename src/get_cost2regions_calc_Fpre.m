function cost = get_cost2regions_calc_Fpre(path, mymodel, model,edata, x, ...
    p_app,gauss_order,prestress_time,eps,changing_matrix,Normalizer,ops_matrix_struct)
    
    global totalRunCount edata_with_Fpre_step ForwardCount
     


    %Creating material parameters vector to sweep on the PVW calculation
    [ground_truth_mat,matparam_sweep,matparam_complete] = SweepMatrix(model,...
        changing_matrix,x,Normalizer,ops_matrix_struct);

    
    % Try-catch block: handle potential failures in FEBio or cost calculation routines.
    try
        % Accumulate prestress calculation for given material parameters updated at every set of parameter.
        % Only update Fpre every 10 evaluations of the cost function
        if ForwardCount==1
            % Always run on the first call, and then every 10th call
            mydir_data = path.data;
            edata = accumulate_Fpre_from_edata(mydir_data, mymodel, ...
                gauss_order, prestress_time, matparam_complete, edata,model);
            edata_with_Fpre_step = edata; % Save the new edata
            
            
            % Uncomment for debug:
            fprintf('[Fpre updated at call %d]\n', totalRunCount);
            fprintf('[Number of Forward count call %d]\n', ForwardCount);
            
        else
            % Use the previous Fpre to save on calculation
            edata = edata_with_Fpre_step;
            % Uncomment for debug:
            fprintf('[Fpre reused at call %d]\n', totalRunCount);
        end


        % Compute the cost function using virtual work integrals.
        [~, ~, ~, cost_func] = calc_virtual_work_variation_integration2(path, mymodel, model, ...
            edata, matparam_complete, matparam_sweep ,ground_truth_mat,p_app,...
            gauss_order,eps,changing_matrix,ops_matrix_struct);

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