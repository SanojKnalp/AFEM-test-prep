classdef ReissnerMindlinSolver < handle
    % REISSNERMINDLINSOLVER Assembles and solves a 2D Reissner-Mindlin Plate.
    
    properties
        mesh            % FEMesh object
        fe              % FEValues object
        E, nu, t, alpha % Material and thickness properties
        q_load          % [3x1] load vector [p_bar; m_x_bar; m_y_bar]
        
        K_global        % Global stiffness matrix
        F_global        % Global load vector
        n_dofs          % Total degrees of freedom
    end
    
    methods
        function obj = ReissnerMindlinSolver(mesh, fe, E, nu, t, alpha, q_load)
            obj.mesh = mesh;
            obj.fe = fe;
            obj.E = E;
            obj.nu = nu;
            obj.t = t;
            obj.alpha = alpha;
            obj.q_load = q_load;
            
            % 3 DOFs per node: [w, beta_x, beta_y]
            obj.n_dofs = size(mesh.nodes, 1) * 3; 
        end
        
        function C = get_C_matrix(obj)
            % GET_C_MATRIX Returns the 5x5 constitutive matrix
            % Rows 1-2: Shear terms
            % Rows 3-5: Bending terms
            
            % Calculate derived properties
            G = obj.E / (2 * (1 + obj.nu));
            D = (obj.E * obj.t^3) / (12 * (1 - obj.nu^2));
            
            C = zeros(5, 5);
            
            % --- YOUR CODE HERE ---
            % Fill in the C matrix based on the provided material law.
            C(1,1) = obj.alpha*G*obj.t; 
            C(2,2) = C(1,1);
            C(3,3) = D;
            C(4,4) = D;
            C(3,4) = obj.nu*D;
            C(4,3) = obj.nu*D;
            C(5,5) = 0.5*D*(1-obj.nu);
        end
        
        function N_mat = get_N_matrix(obj, N_shape)
            % GET_N_MATRIX Assembles the 3x(3n) shape function matrix (N_f)
            
            n_nodes = length(N_shape);
            N_mat = zeros(3, 3 * n_nodes);
            
            % --- YOUR CODE HERE ---
            % Interleave N_shape into the 3x3n matrix.
            for i = 1:1:n_nodes
                N_mat(1,3*i-2) = N_shape(i);
                N_mat(2,3*i-1) = N_shape(i);
                N_mat(3, 3*i) = N_shape(i);
            end
            
        end
        
        function B = get_B_matrix(obj, N_shape, dNdx, dNdy)
            % GET_B_MATRIX Assembles the 5x(3n) kinematic operator matrix
            
            n_nodes = length(N_shape);
            B = zeros(5, 3 * n_nodes);
            
            % --- YOUR CODE HERE ---
            % Assemble the nodal B-blocks. For each node i, insert the 5x3 block B_i.
            for i = 1:1:n_nodes
                B(1, 3*i-2) = dNdx(i);
                B(1, 3*i) = N_shape(i);
                B(2, 3*i-2) = dNdy(i);
                B(2, 3*i-1) = -N_shape(i);
                B(3, 3*i) = dNdx(i);
                B(4, 3*i-1) = -dNdy(i);
                B(5, 3*i-1) = -dNdx(i);
                B(5, 3*i) = dNdy(i);
            end
            
        end
        
        function assemble(obj)
            % ASSEMBLE Allocates and builds the global system.
            obj.K_global = zeros(obj.n_dofs, obj.n_dofs);
            obj.F_global = zeros(obj.n_dofs, 1);
            
            n_elements = size(obj.mesh.elements, 1);
            n_nodes_per_elem = obj.fe.n_nodes;
            
            % 2D Quadrature (Tensor product of 1D rules)
            [xi_1d, w_1d] = quadrature_1d(obj.fe.degree + 1);
            [XI, ETA] = ndgrid(xi_1d, xi_1d);
            [WX, WY] = ndgrid(w_1d, w_1d);
            xi_q = XI(:); eta_q = ETA(:);
            w_q = WX(:) .* WY(:);
            
            C = obj.get_C_matrix();
            
            for e = 1:n_elements
                obj.fe.reinit(e);
                local_nodes = obj.mesh.get_local_nodes(e);
                
                % Map local nodes to global DOFs [w_1, bx_1, by_1, w_2, bx_2, by_2, ...]
                dof_indices = zeros(1, 3 * n_nodes_per_elem);
                for i = 1:n_nodes_per_elem
                    dof_indices(3*i - 2) = 3 * local_nodes(i) - 2; % w
                    dof_indices(3*i - 1) = 3 * local_nodes(i) - 1; % beta_x
                    dof_indices(3*i)     = 3 * local_nodes(i);     % beta_y
                end
                
                K_e = zeros(3 * n_nodes_per_elem, 3 * n_nodes_per_elem);
                F_e = zeros(3 * n_nodes_per_elem, 1);
                
                % Integration loop
                for q = 1:length(w_q)
                    param_coord = [xi_q(q), eta_q(q)];
                    weight = w_q(q);
                    
                    JxW   = obj.fe.get_JxW(param_coord, weight);
                    N     = obj.fe.get_shape(param_coord);
                    dNdx  = obj.fe.get_dNdx(param_coord);
                    dNdy  = obj.fe.get_dNdy(param_coord);
                    
                    B     = obj.get_B_matrix(N, dNdx, dNdy);
                    N_mat = obj.get_N_matrix(N);
                    
                    % Stiffness and load integration
                    K_e = K_e + (B' * C * B) * JxW;
                    F_e = F_e + (N_mat' * obj.q_load) * JxW;
                end
                
                % Global assembly
                obj.K_global(dof_indices, dof_indices) = obj.K_global(dof_indices, dof_indices) + K_e;
                obj.F_global(dof_indices) = obj.F_global(dof_indices) + F_e;
            end
        end
    end
end