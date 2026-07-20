% TEST_EX07_BC_AND_SHEAR
clear; clc;
addpath('../exercises');

fprintf('Running Updated Test Suite (BC Solvers & Shear Forces)...\n\n');
tol_zero = 1e-5;

%% 1. Setup Physical Domain & Mesh
E = 1e6; I = 1/12; G = 5e5; A = 1.0; k_s = 5/6; 
q_load = [-10; 0]; 

% Generate a 5-element mesh over a beam of length 1.0
mesh = FEMesh.generate(1, 5, 1); 
fe = FEValues(mesh, 1, 1);

% Initialize and assemble global arrays
solver = TimoshenkoSolver(mesh, fe, E, I, G, A, k_s, q_load);
solver.assemble();

%% 2. Test Boundary Condition Solvers
% Execute Penalty Method
solver.solve_penalty();
u_penalty = solver.u_global;

% Execute Lagrange Multiplier Method
solver.solve_lagrange_multiplier();
u_lagrange = solver.u_global;

% Execute Direct Elimination Method
solver.solve_direct();
u_direct = solver.u_global;

% Assertions: Verify left-hand clamped BCs (w_1 = 0, beta_1 = 0)
assert(abs(u_penalty(1)) < tol_zero && abs(u_penalty(2)) < tol_zero, 'Penalty: BCs not enforced at left end!');
assert(abs(u_lagrange(1)) < tol_zero && abs(u_lagrange(2)) < tol_zero, 'Lagrange: BCs not enforced at left end!');
assert(abs(u_direct(1)) < tol_zero && abs(u_direct(2)) < tol_zero, 'Direct: BCs not enforced at left end!');

% Assertions: Compare methods against each other
diff_lagrange = max(abs(u_direct - u_lagrange));
assert(diff_lagrange < 1e-12, 'Mismatch between Direct and Lagrange Multiplier methods!');

diff_penalty = max(abs(u_direct - u_penalty));
assert(diff_penalty < 1e-4, 'Mismatch between Direct and Penalty methods (error too high)!');

disp('✅ Boundary condition tests passed! All three methods match.');

%% 3. Test Shear Force Evaluation & Plotting
disp('Evaluating and plotting shear forces...');
num_points = 20;

% Call the updated method using the direct solver results
Q_pts = solver.evaluate_and_plot_shear_forces(num_points);

% Assertions: Verify output integrity
assert(length(Q_pts) == num_points, 'Shear force vector has incorrect length!');
assert(~any(isnan(Q_pts)), 'Shear force vector contains NaN values!');

disp('✅ Shear force evaluation passed! Check the generated figure.');
fprintf('\nAll active tests completed successfully!\n');