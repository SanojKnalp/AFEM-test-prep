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

d

    if degree == 1
        % --- YOUR CODE HERE FOR BILINEAR (4 nodes) ---
        N      = zeros(1, 4);
        dNdxi  = zeros(1, 4);
        dNdeta = zeros(1, 4);
        
    elseif degree == 2
        % --- YOUR CODE HERE FOR BIQUADRATIC (9 nodes) ---
        N      = zeros(1, 9);
        dNdxi  = zeros(1, 9);
        dNdeta = zeros(1, 9);
        
    else
        error('Only degree 1 (bilinear) and 2 (biquadratic) are supported.');
    end
end