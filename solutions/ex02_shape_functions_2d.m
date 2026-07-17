function [N, dNdxi, dNdeta] = ex02_shape_functions_2d(degree, xi, eta)
% EX02_SHAPE_FUNCTIONS_2D Evaluates 2D tensor-product shape functions.
%
% Inputs:
%   degree : integer (1 for bilinear 4-node, 2 for biquadratic 9-node)
%   xi     : double (parametric coordinate in [-1, 1])
%   eta    : double (parametric coordinate in [-1, 1])
%
% Outputs:
%   N      : row vector containing shape function values
%   dNdxi  : row vector containing derivatives w.r.t. xi
%   dNdeta : row vector containing derivatives w.r.t. eta

[Nxi, dNxi] = ex01_shape_functions(degree, xi);
[Neta, dNeta] = ex01_shape_functions(degree, eta);

    if degree == 1
        % --- YOUR CODE HERE FOR BILINEAR (4 nodes) ---
        N      = [Nxi(1)*Neta(1), Nxi(2)*Neta(1), Nxi(1)*Neta(2), Nxi(2)*Neta(2)];
        dNdxi  = [dNxi(1)*Neta(1), dNxi(2)*Neta(1), dNxi(1)*Neta(2), dNxi(2)*Neta(2)];
        dNdeta = [Nxi(1)*dNeta(1), Nxi(2)*dNeta(1), Nxi(1)*dNeta(2), Nxi(2)*dNeta(2)];
        
    elseif degree == 2
        % --- YOUR CODE HERE FOR BIQUADRATIC (9 nodes) ---
        N      = [Nxi(1)*Neta(1), Nxi(2)*Neta(1), Nxi(3)*Neta(1), Nxi(1)*Neta(2), Nxi(2)*Neta(2), Nxi(3)*Neta(2), Nxi(1)*Neta(3), Nxi(2)*Neta(3), Nxi(3)*Neta(3)];
        dNdxi  = [dNxi(1)*Neta(1), dNxi(2)*Neta(1), dNxi(3)*Neta(1), dNxi(1)*Neta(2), dNxi(2)*Neta(2), dNxi(3)*Neta(2), dNxi(1)*Neta(3), dNxi(2)*Neta(3), dNxi(3)*Neta(3)];
        dNdeta = [Nxi(1)*dNeta(1), Nxi(2)*dNeta(1), Nxi(3)*dNeta(1), Nxi(1)*dNeta(2), Nxi(2)*dNeta(2), Nxi(3)*dNeta(2), Nxi(1)*dNeta(3), Nxi(2)*dNeta(3), Nxi(3)*dNeta(3)];
        
    else
        error('Only degree 1 (bilinear) and 2 (biquadratic) are supported.');
    end
end