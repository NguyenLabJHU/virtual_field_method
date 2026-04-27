function [N, dN] = evaluateShapeFunctions(xi_vec)
% EVALUATESHAPEFUNCTIONS - Evaluates shape functions for Quad4 or Hex8
% xi_vec: [xi, eta] for 2D or [xi, eta, zeta] for 3D
% N: Shape function values
% dN: Cell array with derivatives {dN_dxi, dN_deta, dN_dzeta}

dim = length(xi_vec);
xi = xi_vec(1); eta = xi_vec(2);

if dim == 2 % QUAD4 (Surface)
    N = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
    
    dN{1} = 0.25 * [-(1-eta);  (1-eta);  (1+eta); -(1+eta)]; % dN/dxi
    dN{2} = 0.25 * [-(1-xi); -(1+xi);  (1+xi);  (1-xi)];  % dN/deta
    
elseif dim == 3 % HEX8 (Volume)
    zeta = xi_vec(3);
    N = 1/8 * [
        (1-xi)*(1-eta)*(1-zeta); (1+xi)*(1-eta)*(1-zeta);
        (1+xi)*(1+eta)*(1-zeta); (1-xi)*(1+eta)*(1-zeta);
        (1-xi)*(1-eta)*(1+zeta); (1+xi)*(1-eta)*(1+zeta);
        (1+xi)*(1+eta)*(1+zeta); (1-xi)*(1+eta)*(1+zeta)];
        
    dN{1} = 1/8 * [-(1-eta)*(1-zeta); (1-eta)*(1-zeta); (1+eta)*(1-zeta); -(1+eta)*(1-zeta); ...
                   -(1-eta)*(1+zeta); (1-eta)*(1+zeta); (1+eta)*(1+zeta); -(1+eta)*(1+zeta)];
    dN{2} = 1/8 * [-(1-xi)*(1-zeta); -(1+xi)*(1-zeta); (1+xi)*(1-zeta); (1-xi)*(1-zeta); ...
                   -(1-xi)*(1+zeta); -(1+xi)*(1+zeta); (1+xi)*(1+zeta); (1-xi)*(1+zeta)];
    dN{3} = 1/8 * [-(1-xi)*(1-eta); -(1+xi)*(1-eta); -(1+xi)*(1+eta); -(1-xi)*(1+eta); ...
                   (1-xi)*(1-eta);  (1+xi)*(1-eta);  (1+xi)*(1+eta);  (1-xi)*(1+eta)];
end
end