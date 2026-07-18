classdef TimoshenkoSolver < handle
    % TIMOSHENKOSOLVER Assembles and solves a 1D Timoshenko Beam.
    
    properties
        mesh            % FEMesh object
        fe              % FEValues object
        E, I, G, A, k_s % Material and cross-section properties
        q_load          % [2x1] load vector [q_bar; m_bar]
        
        K_global        % Global stiffness matrix
        F_global        % Global load vector
        n_dofs          % Total degrees of freedom
    end
    
    methods
        function obj = TimoshenkoSolver(mesh, fe, E, I, G, A, k_s, q_load)
            obj.mesh = mesh;
            obj.fe = fe;
            obj.E = E;
            obj.I = I;
            obj.G = G;
            obj.A = A;
            obj.k_s = k_s; % This is your alpha
            obj.q_load = q_load;
            
            % 2 DOFs per node: [w, beta]
            obj.n_dofs = size(mesh.nodes, 1) * 2; 
        end
        
        function C = get_C_matrix(obj)
            % GET_C_MATRIX Returns the 2x2 constitutive matrix
            % Row 1: Shear stiffness (alpha * G * A)
            % Row 2: Bending stiffness (EI)
            
            % --- YOUR CODE HERE ---
            C = zeros(2, 2);
            C(1,1) = obj.k_s*obj.G*obj.A;
            C(2,2) = obj.E*obj.I;
        end
        
        function N_mat = get_N_matrix(obj, N_shape)
            % GET_N_MATRIX Assembles the 2x(2n) shape function matrix (N_u)
            % Row 1: Transverse load shape functions
            % Row 2: Moment load shape functions
            
            n_nodes = length(N_shape);
            N_mat = zeros(2, 2 * n_nodes);
            
            % --- YOUR CODE HERE ---
            for i = 1:1: n_nodes
                N_mat(1, 2*i-1) = N_shape(i);
                N_mat(2, 2*i) = N_shape(i);
            end
            
        end
        
        function B = get_B_matrix(obj, N_shape, dNdx)
            % GET_B_MATRIX Assembles the kinematic operator matrix.
            % B should be a [2 x (2*n_nodes)] matrix.
            % Row 1: Shear operator (gamma = w_{,x} + beta)
            % Row 2: Bending operator (kappa = beta_{,x})
            
            n_nodes = length(N_shape);
            B = zeros(2, 2 * n_nodes);
            
            % --- YOUR CODE HERE ---
            for i = 1:1: n_nodes
                B(1,2*i-1 ) = dNdx(i);
                B(1, 2*i) = N_shape(i);
                B(2, 2*i) = dNdx(i);
            end
            
        end
        
        function assemble(obj)
            % ASSEMBLE Allocates and builds the global system.
            obj.K_global = zeros(obj.n_dofs, obj.n_dofs);
            obj.F_global = zeros(obj.n_dofs, 1);
            
            n_elements = size(obj.mesh.elements, 1);
            n_nodes_per_elem = obj.fe.n_nodes;
            
            [xi_q, w_q] = quadrature_1d(obj.fe.degree + 1);
            C = obj.get_C_matrix();
            
            for e = 1:n_elements
                obj.fe.reinit(e);
                local_nodes = obj.mesh.get_local_nodes(e);
                
                dof_indices = zeros(1, 2 * n_nodes_per_elem);
                for i = 1:n_nodes_per_elem
                    dof_indices(2*i - 1) = 2 * local_nodes(i) - 1; % w
                    dof_indices(2*i)     = 2 * local_nodes(i);     % beta
                end
                
                K_e = zeros(2 * n_nodes_per_elem, 2 * n_nodes_per_elem);
                F_e = zeros(2 * n_nodes_per_elem, 1);
                
                for q = 1:length(w_q)
                    param_coord = xi_q(q);
                    weight = w_q(q);
                    
                    JxW   = obj.fe.get_JxW(param_coord, weight);
                    N     = obj.fe.get_shape(param_coord);
                    dNdx  = obj.fe.get_dNdx(param_coord);
                    
                    B     = obj.get_B_matrix(N, dNdx);
                    N_mat = obj.get_N_matrix(N);
                    
                    K_e = K_e + (B' * C * B) * JxW;
                    F_e = F_e + (N_mat' * obj.q_load) * JxW;
                end
                
                obj.K_global(dof_indices, dof_indices) = obj.K_global(dof_indices, dof_indices) + K_e;
                obj.F_global(dof_indices) = obj.F_global(dof_indices) + F_e;
            end
        end
    end
end