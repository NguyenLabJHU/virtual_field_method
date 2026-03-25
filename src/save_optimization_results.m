function save_optimization_results(dataPath, t_elapsed, start_points, fvals, x_opts)
% SAVE_OPTIMIZATION_RESULTS Exports VFM optimization data to a formatted text file.
%
% Inputs:
%   dataPath     - Path to the 'data' folder
%   t_elapsed    - Total convergence time in seconds
%   start_points - Matrix of initial guesses
%   fvals        - Vector of final cost function values
%   x_opts       - Matrix of optimized parameters

    % Define the output file name
    outputFileName = fullfile(dataPath, 'Optimization_Log.txt');

    % Open the file for writing ('w' overwrites existing content)
    fileID = fopen(outputFileName, 'w');

    if fileID == -1
        error('Could not open file for writing. Check folder permissions.');
    end

    % --- Write Header Section ---
    fprintf(fileID, '====================================================\n');
    fprintf(fileID, 'VFM OPTIMIZATION RESULTS - OCULAR BIOMECHANICS\n');
    fprintf(fileID, 'Date: %s\n', datetime('now'));
    fprintf(fileID, 'Total Convergence Time: %d seconds\n', t_elapsed);
    fprintf(fileID, '====================================================\n\n');

    % --- Write Detailed Table ---
    fprintf(fileID, 'Detailed Start Point Analysis:\n');
    fprintf(fileID, '%-8s | %-45s | %-12s | %-45s\n', 'Start #', 'Initial x0', 'Cost (fval)', 'Optimal x');
    fprintf(fileID, '%s\n', repmat('-', 1, 120));

    [n_start, ~] = size(start_points);

    for i = 1:n_start
        % Format arrays as strings for a cleaner table look
        x0_str = sprintf('[%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]', start_points(i,:));
        xopt_str = sprintf('[%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]', x_opts(i,:));
        
        fprintf(fileID, '%-8d | %-45s | %-12.4e | %-45s\n', ...
            i, x0_str, fvals(i), xopt_str);
    end

    % --- Append Raw Matrix ---
    fprintf(fileID, '\n\n--- Raw Optimal Matrix (x_opts) ---\n');
    fclose(fileID);

    % Use writematrix to append the numeric data at the end
    writematrix(x_opts, outputFileName, 'WriteMode', 'append', 'Delimiter', 'tab');

    fprintf('Optimization log successfully saved to: %s\n', outputFileName);
end