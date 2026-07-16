function [xi, w] = quadrature_1d(n_points)
% QUADRATURE_1D Returns Gauss-Legendre points and weights on reference interval [-1, 1]
%
% Inputs:
%   n_points : number of Gauss points (1, 2, or 3)
%
% Outputs:
%   xi       : column vector of coordinates
%   w        : column vector of weights

    if n_points == 1
        xi = 0.0;
        w  = 2.0;
    elseif n_points == 2
        xi = [-1/sqrt(3); 1/sqrt(3)];
        w  = [1.0; 1.0];
    elseif n_points == 3
        xi = [-sqrt(3/5); 0.0; sqrt(3/5)];
        w  = [5/9; 8/9; 5/9];
    else
        error('Unsupported number of quadrature points.');
    end
end