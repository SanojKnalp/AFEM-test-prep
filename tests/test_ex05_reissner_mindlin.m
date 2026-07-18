% TEST_EX06_REISSNER_MINDLIN
clear; clc;
addpath('../exercises');

fprintf('Running Test Suite for Exercise 6 (RM Plate Assembly)...\n\n');
tol = 1e-10;

% Physical parameters
E = 1e6; nu = 0.3; t = 0.1; alpha = 5/6;
q_load = [-10; 0; 0]; % [p_bar; m_x_bar; m_y_bar]

% Derived parameters for testing
G = E / (2 * (1 + nu));
D = (E * t^3) / (12 * (1 - nu^2));

% Create a single-element 2D mesh (Area = 1.0)
mesh = FEMesh.generate(2, 1, 1);
fe = FEValues(mesh, 2, 1);

% Initialize solver
solver = ReissnerMindlinSolver(mesh, fe, E, nu, t, alpha, q_load);

%% Test 1: C-Matrix
C = solver.get_C_matrix();
expected_C = [alpha*G*t, 0, 0, 0, 0;
              0, alpha*G*t, 0, 0, 0;
              0, 0, D, nu*D, 0;
              0, 0, nu*D, D, 0;
              0, 0, 0, 0, 0.5*D*(1 - nu)];
              
assert(all(abs(C(:) - expected_C(:)) < tol), 'C-Matrix is incorrect!');

%% Test 2: N-Matrix and B-Matrix logic
% Check at center point (xi = 0, eta = 0)
solver.fe.reinit(1);
param_coords = [0.0, 0.0];
N = solver.fe.get_shape(param_coords);
dNdx = solver.fe.get_dNdx(param_coords); % Area = 1.0 -> J = diag(0.5, 0.5)
dNdy = solver.fe.get_dNdy(param_coords);

% Test N-Matrix Structure
N_mat = solver.get_N_matrix(N);
expected_N_mat = zeros(3, 12);
for i = 1:4
    expected_N_mat(1, 3*i - 2) = N(i);
    expected_N_mat(2, 3*i - 1) = N(i);
    expected_N_mat(3, 3*i)     = N(i);
end
assert(all(abs(N_mat(:) - expected_N_mat(:)) < tol), 'N-Matrix mapping is incorrect!');

% Test B-Matrix Structure (Checking just Node 1 block for simplicity)
B = solver.get_B_matrix(N, dNdx, dNdy);
B_node1 = B(:, 1:3);
expected_B_node1 = [dNdx(1),       0,    N(1);
                    dNdy(1),   -N(1),       0;
                          0,       0, dNdx(1);
                          0, -dNdy(1),       0;
                          0, -dNdx(1), dNdy(1)];

assert(all(abs(B_node1(:) - expected_B_node1(:)) < tol), 'B-Matrix structure is incorrect!');

fprintf('All RM Plate matrix tests passed! Ready for full plate assembly.\n');