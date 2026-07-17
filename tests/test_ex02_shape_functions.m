% TEST_EX02_SHAPE_FUNCTIONS_2D
clear; clc;
addpath('../exercises');

fprintf('Running Test Suite for Exercise 2...\n\n');
tol = 1e-12;

%% Test 1: Bilinear Element (Q1)
xi_test = 0.2; eta_test = -0.4;
[N_q1, dNdxi_q1, dNdeta_q1] = ex02_shape_functions_2d(1, xi_test, eta_test);

% Partition of unity
assert(abs(sum(N_q1) - 1.0) < tol, 'Q1: Partition of unity failed!');
assert(abs(sum(dNdxi_q1)) < tol, 'Q1: Sum of dNdxi must be 0!');
assert(abs(sum(dNdeta_q1)) < tol, 'Q1: Sum of dNdeta must be 0!');

% Kronecker Delta at Node 3 (-1, 1)
[N_n3, ~, ~] = ex02_shape_functions_2d(1, -1.0, 1.0);
expected_N3 = [0, 0, 1, 0];
assert(all(abs(N_n3 - expected_N3) < tol), 'Q1: Kronecker Delta at Node 3 failed!');

%% Test 2: Biquadratic Element (Q2)
xi_test2 = 0.5; eta_test2 = 0.5;
[N_q2, dNdxi_q2, dNdeta_q2] = ex02_shape_functions_2d(2, xi_test2, eta_test2);

% Partition of unity
assert(abs(sum(N_q2) - 1.0) < tol, 'Q2: Partition of unity failed!');
assert(abs(sum(dNdxi_q2)) < tol, 'Q2: Sum of dNdxi must be 0!');
assert(abs(sum(dNdeta_q2)) < tol, 'Q2: Sum of dNdeta must be 0!');

% Kronecker Delta at Node 5 (0, 0)
[N_n5, dNdxi_n5, dNdeta_n5] = ex02_shape_functions_2d(2, 0.0, 0.0);
expected_N5 = [0, 0, 0, 0, 1, 0, 0, 0, 0];
assert(all(abs(N_n5 - expected_N5) < tol), 'Q2: Kronecker Delta at Node 5 failed!');

% Derivative check at the center (0,0) - Node 5 is a bubble function (1-xi^2)*(1-eta^2)
% Its derivative w.r.t xi is -2*xi*(1-eta^2), which is 0 at the origin.
assert(abs(dNdxi_n5(5)) < tol, 'Q2: Derivative at center is incorrect!');

fprintf('All tests passed successfully! 2D Base is established.\n');