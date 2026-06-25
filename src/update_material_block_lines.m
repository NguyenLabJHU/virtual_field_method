function lines = update_material_block_lines(lines, matparam_complete, model)
    nummat = numel(model.material_model);
    for matn = 1:nummat
        mat_info = model.material_model(matn);
        mat_id_pat = ['<material id="' num2str(matn) '"'];
        
        for i = 1:length(lines)
            if contains(lines{i}, mat_id_pat)
                block_start = i; 
                % mat_type here can be 'prestress material' or 'coupled Mooney-Rivlin'
                % depending on how you populated model.material_model
                mat_type = lower(mat_info.model); 
                
                param_fields = fieldnames(mat_info.parameters);
                for p = 1:numel(param_fields)
                    param_name = param_fields{p};
                    prop = mat_change2prop(mat_type, lower(param_name));
                    
                    if isempty(prop), continue; end
                    
                    % We look for the tag (e.g., <c1>, <k>, <E>)
                    tag = (param_name);
                    tag_pat = ['<' tag '>'];
                    
                    for k = block_start+1:length(lines)
                        % Stops if it reaches the end of the material block (safeguard)
                        if contains((lines{k}), '</material>'), break; end
                        
                        if contains((lower(lines{k})), lower(tag_pat))
                            % EXTRA: Preserve the original line indentation
                            spaces = regexp(lines{k}, '^\s*', 'match');
                            if isempty(spaces), spaces = {''}; end
                            
                            % Updates while keeping the original format
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