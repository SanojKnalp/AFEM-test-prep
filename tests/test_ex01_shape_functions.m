% TEST_EX01_SHAPE_FUNCTIONS
% Run this script to verify your implementation in exercises/ex01_shape_functions.m
clear; clc;

% Add exercises folder to path so MATLAB can find your function
addpath('../solutions');

fprintf('Running Test Suite for Exercise 1...\n\n');

%% Test 1: Linear Element Partition of Unity & Kronecker Delta
tol = 1e-12;
xi_test = -0.5;

[N_lin, dNdxi_lin] = ex01_shape_functions(1, xi_test);

% Partition of unity: sum(N) must equal 1
assert(abs(sum(N_lin) - 1.0) < tol, 'Linear: Partition of unity failed!');
% Sum of derivatives must equal 0
assert(abs(sum(dNdxi_lin)) < tol, 'Linear: Sum of derivatives must be 0!');

% Kronecker Delta property at nodes
[N_node1, ~] = ex01_shape_functions(1, -1.0);
[N_node2, ~] = ex01_shape_functions(1,  1.0);
assert(abs(N_node1(1) - 1.0) < tol && abs(N_node1(2) - 0.0) < tol, 'Linear Node 1 Kronecker Delta failed!');
assert(abs(N_node2(1) - 0.0) < tol && abs(N_node2(2) - 1.0) < tol, 'Linear Node 2 Kronecker Delta failed!');


%% Test 2: Quadratic Element Partition of Unity & Kronecker Delta
xi_test_quad = 0.25;
[N_quad, dNdxi_quad] = ex01_shape_functions(2, xi_test_quad);

% Partition of unity
assert(abs(sum(N_quad) - 1.0) < tol, 'Quadratic: Partition of unity failed!');
assert(abs(sum(dNdxi_quad)) < tol, 'Quadratic: Sum of derivatives must be 0!');

% Kronecker Delta property at nodes
[N_q1, ~] = ex01_shape_functions(2, -1.0);
[N_q2, ~] = ex01_shape_functions(2,  0.0);
[N_q3, ~] = ex01_shape_functions(2,  1.0);

assert(all(abs(N_q1 - [1, 0, 0]) < tol), 'Quadratic Node 1 Kronecker Delta failed!');
assert(all(abs(N_q2 - [0, 1, 0]) < tol), 'Quadratic Node 2 Kronecker Delta failed!');
assert(all(abs(N_q3 - [0, 0, 1]) < tol), 'Quadratic Node 3 Kronecker Delta failed!');


%% Test 3: Specific numerical value check
% Check quadratic derivatives at midpoint xi = 0
[~, dNdxi_mid] = ex01_shape_functions(2, 0.0);
% For N1 = 0.5*xi*(xi-1), dN1/dxi at xi=0 is -0.5
% For N2 = 1 - xi^2, dN2/dxi at xi=0 is 0
% For N3 = 0.5*xi*(xi+1), dN3/dxi at xi=0 is 0.5
expected_derivatives = [-0.5, 0.0, 0.5];
assert(all(abs(dNdxi_mid - expected_derivatives) < tol), 'Quadratic derivatives at xi=0 are incorrect!');

fprintf('All tests passed successfully! Perfect job.\n');