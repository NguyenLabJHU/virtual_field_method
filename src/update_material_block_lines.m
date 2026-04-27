% function lines = update_material_block_lines(lines, matparam_complete, model)
% % Generic routine to update XML lines for material blocks based on parameter struct
% 
% nummat = numel(model.material_model);
% 
% for matn = 1:nummat
%     mat_info = model.material_model(matn);
% 
%     % Find the block for this material ID
%     mat_id_pat = ['<material id="' num2str(matn) '"'];
%     for i = 1:length(lines)
%         if contains(lines{i}, mat_id_pat)
%             % Find elastic block and get type
%             elastic_start = [];
%             mat_type = lower(mat_info.model);
%             % Find start of elastic block for this material
%             for j = i+1:length(lines)
%                 if contains(lower(lines{j}), '<elastic type=')
%                     elastic_start = j;
%                     break;
%                 end
%             end
%             if isempty(elastic_start)
%                 continue
%             end
% 
%             % Now search all fields of this material and update accordingly
%             param_struct = mat_info.parameters;
%             param_fields = fieldnames(param_struct);
%             for p = 1:numel(param_fields)
%                 param_name = param_fields{p};           % E.g., 'C1', 'K', etc
%                 param_val  = param_struct.(param_name); % numeric value
% 
%                 % Use your previous lookup function to find matrix column
%                 prop = mat_change2prop(mat_type, lower(param_name)); % returns column
%                 if isempty(prop)
%                     %warning('Parameter %s in material %s not recognized!', param_name, mat_type);
%                     continue;
%                 end
% 
%                 % Search lines in the XML block for this parameter tag
%                 tag = lower(param_name);
%                 tag_pat = ['<' tag '>'];
%                 for k = elastic_start+1:min(elastic_start+15, length(lines)) % scan a bit ahead in the block
%                     if contains(lower(lines{k}), tag_pat)
%                         lines{k} = sprintf('<%s>%g</%s>', tag, matparam_complete(matn,prop), tag);
%                         break
%                     end
%                 end
%             end % param loop
% 
%             break % Only one block per material
%         end
%     end
% end
% end


function lines = update_material_block_lines(lines, matparam_complete, model)
    nummat = numel(model.material_model);
    for matn = 1:nummat
        mat_info = model.material_model(matn);
        mat_id_pat = ['<material id="' num2str(matn) '"'];
        
        for i = 1:length(lines)
            if contains(lines{i}, mat_id_pat)
                block_start = i; 
                % mat_type aqui pode ser 'prestress material' ou 'coupled Mooney-Rivlin'
                % dependendo de como você preencheu o model.material_model
                mat_type = lower(mat_info.model); 
                
                param_fields = fieldnames(mat_info.parameters);
                for p = 1:numel(param_fields)
                    param_name = param_fields{p};
                    prop = mat_change2prop(mat_type, lower(param_name));
                    
                    if isempty(prop), continue; end
                    
                    % Procuramos a tag (ex: <c1>, <k>, <E>)
                    tag = (param_name);
                    tag_pat = ['<' tag '>'];
                    
                    for k = block_start+1:length(lines)
                        % Para se chegar no fim do material (proteção)
                        if contains((lines{k}), '</material>'), break; end
                        
                        if contains((lower(lines{k})), lower(tag_pat))
                            % EXTRA: Preservar a indentação original da linha
                            spaces = regexp(lines{k}, '^\s*', 'match');
                            if isempty(spaces), spaces = {''}; end
                            
                            % Atualiza mantendo o formato original
                            expr = '(\s*)<([^>]+)>(.*?)</([^>]+)>';
                            lines{k} = regexprep( lines{k}, expr, ...
                                ['$1<$2>' num2str(matparam_complete(matn, prop)) '</$4>'] );
                            break
                        end
                    end
                end
                break 
            end
        end
    end
end