function febio_bin = get_febio_path()
% GET_FEBIO_PATH Automatically locates the FEBio 4 executable.


    % 1. Check for Environment Variable first (Good practice for clusters)
    febio_bin = getenv('FEBIO_PATH');
    if ~isempty(febio_bin) && isfile(febio_bin)
        return;
    end

    % 2. System Path Search
    if ispc % Windows
        [status, cmdOut] = system('where febio4');
    else % Linux / Mac
        [status, cmdOut] = system('which febio4');
    end

    if status == 0
        paths = strsplit(strtrim(cmdOut), char(10)); 
        febio_bin = paths{1}; 
        return;
    end

    % 3. Fallback: Specific paths for your current machines
    if ispc
        % Your specific Windows path in Dell
        febio_bin = 'C:\Program Files\FEBioStudio\bin\febio4.exe';
    else
        % Your specific path on the HPC server
        febio_bin = '/home/msampai4/FEBio/febio_src/build/bin/febio4';
    end

    % Final check to ensure the file actually exists
    if ~isfile(febio_bin)
        warning('FEBio executable not found at: %s. Please check your installation path.', febio_bin);
    end
end