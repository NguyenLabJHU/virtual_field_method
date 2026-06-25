function edata_noisy = pure_gaussian_noise_edata(edata, model, noise_percent, sigma_additive)
    % pure_gaussian_noise_checked: Aplica ruído branco gaussiano puro
    % e verifica o Jacobiano para garantir que nenhum elemento inverta.
    
    edata_noisy = edata;
    nsteps = size(edata.steps, 2);
    coords_ref = edata.steps{1, 1}.results.ecoords;
    [num_nodes, dim] = size(coords_ref);
    
    % Extrai a conectividade (assumindo elementos C3D8 - 8 nós)
    conn = model.elements(:, 2:9);

    % --- 1. Definindo a Escala Baseada no Deslocamento ---
    all_disps = [];
    for s = 1:nsteps
        all_disps = [all_disps; edata.steps{s}.results.edisp]; 
    end
    max_u = max(abs(all_disps(:)));
    noise_std = (noise_percent/100) * max_u + sigma_additive;

    % --- 2. Injeção de Ruído com Válvula de Segurança ---
    fprintf('Iniciando injeção de ruído Gaussiano puro com checagem Jacobiana...\n');
    for step_idx = 1:nsteps
        ecoords_orig = edata.steps{1, step_idx}.results.ecoords;
        
        success = false;
        attempts = 0;
        current_scale = 1.0; % Começa com 100% da escala de ruído solicitada

        % Tenta até 10 vezes reduzir o ruído caso os elementos invertam
        while ~success && attempts < 10
            % Gera ruído branco PURO (independente para cada nó)
            noise_raw = (noise_std * current_scale) .* randn(num_nodes, dim);
            
            % Aplica nas coordenadas (SEM SMOOTH)
            test_coords = ecoords_orig + noise_raw;
            
            % --- CHECAGEM DO JACOBIANO ---
            if check_jacobian(test_coords, conn)
                % Passou no teste! Aceita as coordenadas e encerra o loop
                edata_noisy.steps{1, step_idx}.results.ecoords = test_coords;
                success = true;
            else
                % Falhou. Aumenta a tentativa e corta o ruído pela metade
                attempts = attempts + 1;
                current_scale = current_scale * 0.5;
                fprintf('Aviso (Passo %d): Elemento inverteu. Reduzindo escala de ruído para %.4f\n', step_idx, current_scale);
            end
        end
        
        % Se mesmo após 10 tentativas a malha quebrar, mantém o dado original
        if ~success
            warning('Passo %d: Impossível injetar ruído sem inverter a malha. Mantendo coordenadas originais.', step_idx);
            edata_noisy.steps{1, step_idx}.results.ecoords = ecoords_orig;
        end
    end

    % --- 3. Consistência Cinemática Final ---
    % Atualiza os deslocamentos (u = x - X) para bater com as coordenadas ruidosas
    coords1_noisy = edata_noisy.steps{1, 1}.results.ecoords;
    for step_idx = 1:nsteps
        curr_noisy = edata_noisy.steps{1, step_idx}.results.ecoords;
        edata_noisy.steps{1, step_idx}.results.edisp = curr_noisy - coords1_noisy;
    end
    
    fprintf('Injeção de ruído concluída.\n');
end

% =========================================================================
% FUNÇÃO LOCAL: CHECAGEM DO JACOBIANO
% =========================================================================
function elementos_ok = check_jacobian(coords, conn)
    % coords: matriz [nós x 3] com as coordenadas atuais (x, y, z)
    % conn:   matriz [elementos x 8] com as conectividades da malha
    
    num_elements = size(conn, 1);
    elementos_ok = true; 
    
    % Derivadas das funções de forma no centro do elemento hexaédrico (C3D8)
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
        nos_elemento = conn(i, :);
        coords_elemento = coords(nos_elemento, :); 
        
        % Matriz Jacobiana no centro do elemento
        J = dN' * coords_elemento;
        
        % Verifica o determinante
        if det(J) <= 0
            elementos_ok = false;
            return; % Sai no primeiro erro encontrado para otimizar velocidade
        end
    end
end