function sigma = computeStress(F, model,mat_idx, imatprop,ielement,gp)
    % Compute Cauchy stress for material in 'model' at index 'mat_idx'
    % F: 3x3 deformation gradient
    % model: struct from your parser (see above)
    % mat_idx: index of material to use

    % Determine material model string
    mat_type = lower(strtrim(model.matmodel(mat_idx)));

    % Assign (default to zero if missing)
    prop1     = imatprop(1); 
    prop2     = imatprop(2);
    prop3     = imatprop(3);
    prop4     = imatprop(4);
    prop5     = imatprop(5);
    prop6     = imatprop(6);
    prop7     = imatprop(7);


    %% --------- Retrieving fiber and mat-axis ----------
    % Fiber direction, if present, otherwise [1;0;0]
    if isfield(model.material_model, 'fibers') && ~isempty(model.material_model(mat_idx).fibers)
        % Find the row where first column matches ielement AND second column matches gp
        idx = find( ...
            model.material_model(mat_idx).fibers(:,1) == ielement & ...
            model.material_model(mat_idx).fibers(:,2) == gp, ...
            1);
        if ~isempty(idx)
            a0 = model.material_model(mat_idx).fibers(idx, 3:end).'; % COLUMN
        else

            % Display info about default assignment:
            fprintf('Default fiber assigned to element %d (material type: %s)\n', ...
            ielement, mat_type);
            
            a0 = [1; 0; 0]; % Default if not found
        end
    else
        a0 = [1; 0; 0];
    end


    % Local material axes (mat_axis), if present, otherwise identity
    if isfield(model.material_model, 'mat_axis') && ~isempty(model.material_model(mat_idx).mat_axis)
        idx = find( ...
            model.material_model(mat_idx).mat_axis(:,1) == ielement & ...
            model.material_model(mat_idx).mat_axis(:,2) == gp, ...
            1);
        if ~isempty(idx)
            Q_flat = model.material_model(mat_idx).mat_axis(idx, 3:end);
            Q = reshape(Q_flat, 3, 3);  % Convert to 3x3 matrix
        else
            fprintf('Default mat_axis (eye(3)) assigned to element %d (material type: %s)\n', ...
            ielement, mat_type);
            Q = eye(3); % Default if not found
        end
    else
        Q = eye(3);
    end


    %% --------- Choosing the material ----------
    % Build matprop struct as expected by the solver
    switch mat_type
        case "linear-elastic"
            matprop.E = prop1;   % E
            matprop.v = prop3;   % poisson
        case "neo-hookean"
            lambda = prop3*prop1/((1+prop3)*(1-2*prop3));
            mu = 0.5*prop1/(1+prop3);
            matprop.mu = mu;
            matprop.k  = lambda;
        case "coupled-mooney-rivlin"
            matprop.C1 = prop1;
            matprop.C2 = prop2;
            matprop.k  = prop3;
        case "mooney-rivlin"
            matprop.C1 = prop1;
            matprop.C2 = prop2;
            matprop.k  = prop3;    
        case "ti-mooney-rivlin"
            matprop.C1 = prop1;
            matprop.C2 = prop2;
            matprop.K  = prop3;
            matprop.c3 = prop4;
            matprop.c4 = prop5;
            matprop.c5 = prop6;
            matprop.flam = prop7;
            matprop.a0 = a0;
        case 'hgo-unconstrained'
            matprop.c      = prop1;
            matprop.k1     = prop2;
            matprop.k2     = prop4;
            matprop.K  = prop3;
            matprop.kappa  = prop5;
            matprop.gdeg   = prop6;
            matprop.Q = Q;

        otherwise
            error('Unknown material model type "%s"!', mat_type);
    end

    % Now call inner stress calculator for a per-material struct:
    sigma = computeStress_core(F, matprop, mat_type);

end

function sigma = computeStress_core(F, matprop, model)
% Internal: calculate stress for a constructed matprop and a string model name

    Id = eye(3);
    J  = det(F);
    if J <= 0, error('Invalid deformation: det(F)=J <= 0.'); end
    B  = F * F.';        % Left Cauchy-Green tensor
    C  = F.' * F;        % Right Cauchy-Green tensor
    B2 = B * B;
    I1 = trace(B);

    switch lower(strrep(model,'_','-'))

        case {'neo-hookean','nh'}
            mu    = matprop.mu;
            kappa = matprop.k;
            sigma = (mu/J) * (B - Id) + (kappa*log(J)/J) * Id;

        case {'coupled-mooney-rivlin','cmr'}
            C1 = matprop.C1;
            C2 = matprop.C2;
            kappa = matprop.k;
            sigma_iso = (2.0/J) * ( B*(C1 + I1*C2) - B2*C2 - Id*(C1 + 2.0*C2) );
            sigma_vol = (kappa * log(J) / J) * Id;        
            sigma = sigma_iso + sigma_vol;

        case {'ti-mooney-rivlin','transversely-isotropic-mooney-rivlin','timr'}
            C1   = matprop.C1;
            C2   = matprop.C2;
            K    = matprop.K;
            a0   = matprop.a0(:);
            c3   = matprop.c3;
            c4   = matprop.c4;
            c5   = matprop.c5;
            flam = matprop.flam;

            n0 = norm(a0);
            if n0 == 0, error('TI-MR: material.a0 must be non-zero.'); end
            a0 = a0 / n0;
            sigma_mat = (2.0/J) * ( B*(C1 + I1*C2) - B2*C2 - Id*(C1 + 2.0*C2) );
            a = F * a0;
            l = norm(a);
            if l > 0, a_hat = a / l; A = a_hat*a_hat.'; else, A = zeros(3); end
            sigma_fiber = zeros(3);
            if l > 1.0
                if l < flam
                    Wl = c3 * (exp(c4*(l - 1.0)) - 1.0);
                else
                    c6 = c3*(exp(c4*(flam - 1.0)) - 1.0) - c5*flam;
                    Wl = c5*l + c6;
                end
                sigma_fiber = (Wl / J) * A;
            end
            sigma_vol = (K * log(J) / J) * Id;
            sigma = sigma_mat + sigma_fiber + sigma_vol;


        case {'holzapfel-unconstrained','hgo-unconstrained'}
            % Parameters
            c     = matprop.c;
            k1    = matprop.k1;
            k2    = matprop.k2;
            kappa = matprop.kappa;
            gdeg  = matprop.gdeg * pi/180; % fiber angle in radians
            K     = matprop.K;
            Q    = matprop.Q;         
            
            % (Optionally get Q, the local coordinate system)
            if isfield(matprop, 'Q')
                n1 = matprop.Q(:,1);
                n2 = matprop.Q(:,2);
            else
                % If only a0 is available, build an orthonormal system
                n1 = a0;
                if abs(a0(1)) < 0.99
                    n2 = cross(a0, [1;0;0]);
                else
                    n2 = cross(a0, [0;1;0]);
                end
                n2 = n2 / norm(n2);
            end
    
            % Build fiber directions rotated by +gdeg and -gdeg
            cg = cos(gdeg); sg = sin(gdeg);
            ar{1} = n1 * cg + n2 * sg;
            ar{2} = n1 * cg - n2 * sg;
    
            J = det(F);
            I1 = trace(C);  % First invariant
    
            % Ground matrix: isotropic Cauchy stress part
            sigma = c * (B - eye(3));
    
            % Volumetric penalty: 0.5 * K * (J^2 - 1) * I
            sigma = sigma + 0.5 * K * (J^2 - 1) * eye(3);
    
            % Fiber contributions: two families rotated within fiber plane
            for i = 1:2
                ari = ar{i};
                ai  = F * ari;
                I4  = ari.' * (C * ari);
                E   = kappa * (I1 - 3) + (1 - 3*kappa)*(I4 - 1);
                if E >= 0
                    h = kappa * B + (1 - 3*kappa) * (ai * ai.');
                    sigma = sigma + h * (2 * k1 * E * exp(k2 * E * E));
                end
            end
    
            % Final Cauchy stress (divide by J)
            sigma = sigma / J;

        case {'linear-elastic','linear elastic','le'}
            E = matprop.E;
            v = matprop.v;

            Ji = 1.0 / J;

            % lame parameters (same FEBio)
            lam = Ji*(v*E/((1+v)*(1-2*v)));
            mu  = Ji*(0.5*E/(1+v));

            % FEBio uses b = LeftCauchyGreen and b2 = b.sqr()
            trE = 0.5*(trace(B) - 3.0);
            sigma = B*(lam*trE - mu) + B2*mu;


        case {'mooney-rivlin','dmr'}
            c1 = matprop.C1;
            c2 = matprop.C2;
            kappa = matprop.k;

            % --- FEBio: DevLeftCauchyGreen ---
            Jm23 = J^(-2/3);
            Bbar  = Jm23 * B;        % deviatoric/isocoric B
            Bbar2 = Bbar * Bbar;

            I1 = trace(Bbar);
            I2 = 0.5*(I1*I1 - trace(Bbar2));

            W1 = c1;
            W2 = c2;

            T = Bbar*(W1 + W2*I1) - Bbar2*W2;

            % deviatoric part: Tdev = T - (tr(T)/3)*I
            Tdev = T - (trace(T)/3)*Id;

            sigma_dev = (2.0/J) * Tdev;

            % FEBio normalmente soma a parte volumétrica em outro lugar
            % Se você quer "total stress" com penalidade volumétrica:
            sigma_vol = (kappa * log(J) / J) * Id;

            sigma = sigma_dev + sigma_vol;    


        otherwise
            error('Unknown material model: %s', model);
    end
end