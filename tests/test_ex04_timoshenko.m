% TEST_EX05_TIMOSHENKO_ASSEMBLY
clear; clc;
addpath('../solutions');

fprintf('Running Test Suite for Exercise 5 (Timoshenko Assembly)...\n\n');
tol = 1e-10;

% Physical parameters
E = 1e6; I = 1/12; G = 5e5; A = 1.0; k_s = 5/6; 
q_load = [-10; 0]; % [q_bar; m_bar]

% Create a single-element mesh (Default length = 1.0)
mesh = FEMesh.generate(1, 1, 1);
fe = FEValues(mesh, 1, 1);

% Initialize solver
solver = TimoshenkoSolver(mesh, fe, E, I, G, A, k_s, q_load);

%% Test 1: C-Matrix
C = solver.get_C_matrix();
expected_C = [k_s*G*A, 0; 0, E*I];
assert(all(abs(C(:) - expected_C(:)) < tol), 'C-Matrix is incorrect!');

%% Test 2: N-Matrix and B-Matrix logic
% Check at xi = 0 (midpoint of the reference element)
solver.fe.reinit(1);
N = solver.fe.get_shape(0);
dNdx = solver.fe.get_dNdx(0); % Length = 1.0 -> J = 0.5 -> dNdx = [-1.0, 1.0]

N_mat = solver.get_N_matrix(N);
% N at xi=0 is [0.5, 0.5]. 
expected_N_mat = [0.5, 0, 0.5, 0;
                  0, 0.5, 0, 0.5];
assert(all(abs(N_mat(:) - expected_N_mat(:)) < tol), 'N-Matrix mapping is incorrect!');

B = solver.get_B_matrix(N, dNdx);
% Row 1 (Shear)   : [dNdx1, N1, dNdx2, N2] = [-1.0, 0.5, 1.0, 0.5]
% Row 2 (Bending) : [0, dNdx1, 0, dNdx2]   = [0, -1.0, 0, 1.0]
expected_B = [-1.0, 0.5, 1.0, 0.5; 
               0, -1.0, 0, 1.0];

assert(all(abs(B(:) - expected_B(:)) < tol), 'B-Matrix is incorrect!');

fprintf('All basic matrix tests passed! Ready for full assembly.\n');