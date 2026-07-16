function [N, dNdxi] = ex01_shape_functions(degree, xi)
% EX01_SHAPE_FUNCTIONS Evaluates 1D shape functions and reference derivatives.
%
% Inputs:
%   degree : integer (1 for linear, 2 for quadratic)
%   xi     : double (parametric coordinate in [-1, 1])
%
% Outputs:
%   N      : row vector [1 x n_nodes] containing shape function values
%   dNdxi  : row vector [1 x n_nodes] containing derivatives w.r.t. xi
%
% Lexicographical node ordering on reference element [-1, 1]:
%   degree = 1:  Node 1 @ xi=-1, Node 2 @ xi=1
%   degree = 2:  Node 1 @ xi=-1, Node 2 @ xi=0, Node 3 @ xi=1

    if degree == 1
        % --- YOUR CODE HERE FOR LINEAR ELEMENT ---
        N     = zeros(1, 2); 
        dNdxi = zeros(1, 2);
        
    elseif degree == 2
        % --- YOUR CODE HERE FOR QUADRATIC ELEMENT ---
        N     = zeros(1, 3);
        dNdxi = zeros(1, 3);
        
    else
        error('Only degree 1 (linear) and degree 2 (quadratic) are supported.');
    end
end