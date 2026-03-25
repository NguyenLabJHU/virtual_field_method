function edata_noisy = dirty_steps_edata(edata, noise_percent, sigma_additive)
% DIRTY_STEPS_EDATA Adds gaussian noise to ecoords for all steps,
% and calculates edisp as the difference between noisy ecoords in each step and step 1.
%
% edata_noisy = dirty_steps_edata(edata, noise_percent, sigma_additive)
% edata: struct with steps{1, step_idx}.results.ecoords (Nx3)
% noise_percent: percentage of noise to add to each step's ecoords
% sigma_additive: minimum additive noise (same units as ecoords)
% Returns: edata_noisy, with noisy ecoords and noisy edisp for all steps

    edata_noisy = edata;
    nsteps = size(edata.steps, 2);

    % Generate noisy ecoords for step 1
    ecoords_step1 = edata.steps{1, 1}.results.ecoords;
    ecoords1_noise = randn(size(ecoords_step1)) .* ...
                     (noise_percent/100 .* abs(ecoords_step1) + sigma_additive);
    ecoords1_noisy = ecoords_step1 + ecoords1_noise;

    % For step 1, displacement is zero (or difference with itself)
    edata_noisy.steps{1, 1}.results.ecoords = ecoords1_noisy;
    edata_noisy.steps{1, 1}.results.edisp = zeros(size(ecoords1_noisy));

    for step_idx = 2:nsteps
        % Get the current ecoords
        ecoords_cur = edata.steps{1, step_idx}.results.ecoords;

        % Add noise to ecoords_cur
        ecoords_cur_noise = randn(size(ecoords_cur)) .* ...
                            (noise_percent/100 .* abs(ecoords_cur) + sigma_additive);
        %ecoords_cur_noisy = ecoords_cur + ecoords_cur_noise;
        ecoords_cur_noisy = ecoords_cur;

        % Compute edisp as di   fference between noisy current and noisy step 1 ecoords
        edisp_cur_noisy = ecoords_cur_noisy - ecoords1_noisy;

        % Store in output struct
        %edata_noisy.steps{1, step_idx}.results.ecoords = ecoords_cur_noisy;
        edata_noisy.steps{1, step_idx}.results.edisp = edisp_cur_noisy;
    end
end