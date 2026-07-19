classdef FEValues < handle
    % FEVALUES On-the-fly FE evaluation class with centralized Jacobian.
    
    properties (SetAccess = private)
        mesh         % Reference to the FEMesh object
        dim          % Dimension (1 or 2)
        degree       % Polynomial degree (1 or 2)
        elem_coords  % [n_nodes x dim] physical coordinates of current element
    end
    
    methods
        function obj = FEValues(mesh, dim, degree)
            obj.mesh = mesh;
            obj.dim = dim;
            obj.degree = degree;
        end
        
        function reinit(obj, elem_id)
            obj.elem_coords = obj.mesh.get_element_coordinates(elem_id);
        end
        
        function N = get_shape(obj, param_coords)
            if obj.dim == 1
                [N, ~] = ex01_shape_functions(obj.degree, param_coords(1));
            elseif obj.dim == 2
                [N, ~, ~] = ex02_shape_functions_2d(obj.degree, param_coords(1), param_coords(2));
            end
        end
        
        function J = get_jacobian(obj, param_coords)
            % GET_JACOBIAN Computes the Jacobian matrix (1x1 for 1D, 2x2 for 2D)
            if obj.dim == 1
                % --- YOUR CODE HERE FOR 1D JACOBIAN ---
                % 1. Get reference derivative dNdxi
                [~,dN] = ex01_shape_functions(obj.degree, param_coords(1));
                % 2. Compute J = dNdxi * elem_coords
                J = 0; 
                for i = 1:1:(obj.degree+1)
                    J = J+ dN(i)*obj.elem_coords(i);
                end
                
            elseif obj.dim == 2
                % --- YOUR CODE HERE FOR 2D JACOBIAN ---
                % 1. Get reference derivatives dNdxi, dNdeta
                [~, dNdxi, dNdeta] = ex02_shape_functions_2d(obj.degree,param_coords(1), param_coords(2));
                % 2. Compute 2x2 matrix J = [dx/dxi, dy/dxi ; dx/deta, dy/deta]
                % Hint: You can do this elegantly via matrix multiplication: 
                % J = [dNdxi; dNdeta] * obj.elem_coords;
                J = [dNdxi; dNdeta] * obj.elem_coords; 

            end
        end
        
        function JxW = get_JxW(obj, param_coords, weight)
            % GET_JXW Computes det(J) * weight
            J = obj.get_jacobian(param_coords);
            
            % --- YOUR CODE HERE ---
            % Compute the determinant of J and multiply by the weight
            JxW = det(J)*weight;
        end
        
        function dNdx = get_dNdx(obj, param_coords)
            % GET_DNDX Computes physical x-derivatives of shape functions
            J = obj.get_jacobian(param_coords);
            
            if obj.dim == 1
                % --- YOUR CODE HERE FOR 1D ---
                [~,dN] = ex01_shape_functions(obj.degree, param_coords(1));
                dNdx = zeros(1, size(obj.elem_coords, 1));
                dNdx(1) = dN(1)/J;
                dNdx(2) = dN(2)/J;
                
            elseif obj.dim == 2
                % --- YOUR CODE HERE FOR 2D ---
                [~, dNdxi, dNdeta] = ex02_shape_functions_2d(obj.degree,param_coords(1), param_coords(2));
                dNdx = zeros(1, size(obj.elem_coords, 1));
                vec = inv(J)*[ dNdxi; dNdeta];
                dNdx = vec(1, :);

            end
        end
        
        function dNdy = get_dNdy(obj, param_coords)
            % GET_DNDY Computes physical y-derivatives of shape functions (2D only)
            if obj.dim == 1
                error('y-derivative is not defined for a 1D problem.');
            elseif obj.dim == 2
                J = obj.get_jacobian(param_coords);
                
                % --- YOUR CODE HERE FOR 2D ---
                [~, dNdxi, dNdeta] = ex02_shape_functions_2d(obj.degree,param_coords(1), param_coords(2));
                dNdy = zeros(1, size(obj.elem_coords, 1));
                vec = inv(J)*[ dNdxi; dNdeta];
                dNdy = vec(2, :);
            end
        end
    end
end