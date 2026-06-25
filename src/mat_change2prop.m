function prop = mat_change2prop(mat_type, mat_change)
%MAT_CHANGE2PROP Get the property index for a given parameter and material type.
%   prop = MAT_CHANGE2PROP(mat_type, mat_change)
%
%   Inputs:
%     mat_type   - String with material type (e.g. 'neo-hookean')
%     mat_change - String with parameter name (e.g. 'c1', 'k', etc)
%
%   Output:
%     prop       - Property index (numeric); [] if not found

    prop = []; % default if not found

    switch lower(mat_type)
        case "neo-hookean"
            switch lower(mat_change)
                case 'mu'
                    prop = 1;
                case 'k'
                    prop = 2;
            end
        case {"mooney-rivlin", "coupled-mooney-rivlin"}
            switch lower(mat_change)
                case 'c1'
                    prop = 1;
                case 'c2'
                    prop = 2;
                case 'k'
                    prop = 3;
            end
        case "ti-mooney-rivlin"
            switch lower(mat_change)
                case 'c1'
                    prop = 1;
                case 'c2'
                    prop = 2;
                case 'k'
                    prop = 3;
                case 'c3'
                    prop = 4;
                case 'c4'
                    prop = 5;
                case 'c5'
                    prop = 6;
                case 'flam'
                    prop = 7;
            end
        case "hgo-unconstrained"
            switch lower(mat_change)
                case 'c'
                    prop = 1;
                case 'k1'
                    prop = 2;
                case 'k'
                    prop = 3;
                case 'k2'
                    prop = 4;
                case 'kappa'
                    prop = 5;
                case 'gamma'
                    prop = 6;
            end
        case {"linear-elastic","linear elastic","le"}
            switch lower(mat_change)
                case 'e'
                    prop = 1;
                case {'v','nu','poisson'}
                    prop = 3;
            end
    end
end