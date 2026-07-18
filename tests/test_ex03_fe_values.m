% TEST_EX04_FEVALUES
clear; clc;
addpath('../exercises');

fprintf('Running Test Suite for Exercise 4 (Dynamic FEValues)...\n\n');
tol = 1e-12;

%% Test 1: 1D Linear Element Evaluation
mesh_1d = FEMesh.generate(1, 1, 1); % Length = 1.0
fe_1d = FEValues(mesh_1d, 1, 1);
fe_1d.reinit(1); 

% Test Jacobian directly
J_1d = fe_1d.get_jacobian(0.0);
assert(abs(J_1d - 0.5) < tol, '1D: get_jacobian is incorrect!');

JxW_1d = fe_1d.get_JxW(0.0, 2.0);
assert(abs(JxW_1d - 1.0) < tol, '1D: get_JxW is incorrect!');

dNdx_1d = fe_1d.get_dNdx(0.5);
assert(all(abs(dNdx_1d - [-1.0, 1.0]) < tol), '1D: get_dNdx is incorrect!');


%% Test 2: 2D Bilinear Element Evaluation
mesh_2d = FEMesh.generate(2, 1, 1); % Area = 1.0
fe_2d = FEValues(mesh_2d, 2, 1);
fe_2d.reinit(1); 

% Test Jacobian directly
J_2d = fe_2d.get_jacobian([0.0, 0.0]);
expected_J2d = [0.5, 0.0; 0.0, 0.5];
assert(all(abs(J_2d(:) - expected_J2d(:)) < tol), '2D: get_jacobian is incorrect!');

JxW_2d = fe_2d.get_JxW([0.0, 0.0], 4.0);
assert(abs(JxW_2d - 1.0) < tol, '2D: get_JxW is incorrect!');

dNdx_2d = fe_2d.get_dNdx([0.0, 0.0]);
dNdy_2d = fe_2d.get_dNdy([0.0, 0.0]);

expected_dNdx = [-0.5, 0.5, -0.5, 0.5];
expected_dNdy = [-0.5, -0.5, 0.5, 0.5];

assert(all(abs(dNdx_2d - expected_dNdx) < tol), '2D: get_dNdx is incorrect!');
assert(all(abs(dNdy_2d - expected_dNdy) < tol), '2D: get_dNdy is incorrect!');

fprintf('All tests passed successfully! Your dynamic FEValues class is fully operational.\n');