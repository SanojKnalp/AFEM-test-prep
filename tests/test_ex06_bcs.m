% TEST_EX07_BC_METHODS
clear; clc;
addpath('../solutions');

fprintf('Running Test Suite for Exercise 7 (Boundary Condition Methods)...\n\n');

% Set up physical parameters for a standard beam
E = 1e6; I = 1/12; G = 5e5; A = 1.0; k_s = 5/6; 
q_load = [-10; 0]; 

% Create a 5-element mesh to ensure we have enough DOFs to test
mesh = FEMesh.generate(1, 5, 1);
fe = FEValues(mesh, 1, 1);

% Initialize and assemble
solver = TimoshenkoSolver(mesh, fe, E, I, G, A, k_s, q_load);
solver.assemble();

%% Solve using all three methods
penalty_val = 1e12;
u_penalty  = solver.solve_penalty(penalty_val);
u_lagrange = solver.solve_lagrange_multiplier();
u_direct   = solver.solve_direct();

%% Testing Asserts

% 1. Check if BCs are actually enforced at Node 1 (w=0, beta=0)
tol_zero = 1e-10;
assert(abs(u_penalty(1)) < tol_zero && abs(u_penalty(2)) < tol_zero, 'Penalty: BCs not enforced at left end!');
assert(abs(u_lagrange(1)) < tol_zero && abs(u_lagrange(2)) < tol_zero, 'Lagrange: BCs not enforced at left end!');
assert(abs(u_direct(1)) < tol_zero && abs(u_direct(2)) < tol_zero, 'Direct: BCs not enforced at left end!');

% 2. Compare Direct and Lagrange (Should be mathematically identical)
diff_lagrange = max(abs(u_direct - u_lagrange));
assert(diff_lagrange < 1e-12, 'Mismatch between Direct and Lagrange Multiplier methods!');

% 3. Compare Direct and Penalty (Will have slight error due to penalty stiffness)
diff_penalty = max(abs(u_direct - u_penalty));
assert(diff_penalty < 1e-4, 'Mismatch between Direct and Penalty methods (error too high)!');

fprintf('All boundary condition tests passed! Deflection fields match perfectly.\n');