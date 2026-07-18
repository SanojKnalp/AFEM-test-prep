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
                % fill 
            elseif obj.dim == 2
               % fill 
            end
        end
        
        function J = get_jacobian(obj, param_coords)
            % GET_JACOBIAN Computes the Jacobian matrix (1x1 for 1D, 2x2 for 2D)
            if obj.dim == 1
                % --- YOUR CODE HERE FOR 1D JACOBIAN ---
                % 1. Get reference derivative dNdxi
                % 2. Compute J = dNdxi * elem_coords
                J = 0; 
                
            elseif obj.dim == 2
                % --- YOUR CODE HERE FOR 2D JACOBIAN ---
                % 1. Get reference derivatives dNdxi, dNdeta
                % 2. Compute 2x2 matrix J = [dx/dxi, dy/dxi ; dx/deta, dy/deta]
                % Hint: You can do this elegantly via matrix multiplication: 
                % J = [dNdxi; dNdeta] * obj.elem_coords;
                J = zeros(2, 2); 
            end
        end
        
        function JxW = get_JxW(obj, param_coords, weight)
            % GET_JXW Computes det(J) * weight
            J = obj.get_jacobian(param_coords);
            
            % --- YOUR CODE HERE ---
            % Compute the determinant of J and multiply by the weight
            JxW = 0; 
        end
        
        function dNdx = get_dNdx(obj, param_coords)
            % GET_DNDX Computes physical x-derivatives of shape functions
            J = obj.get_jacobian(param_coords);
            
            if obj.dim == 1
                % --- YOUR CODE HERE FOR 1D ---
                dNdx = zeros(1, size(obj.elem_coords, 1));
                
            elseif obj.dim == 2
                % --- YOUR CODE HERE FOR 2D ---
                dNdx = zeros(1, size(obj.elem_coords, 1));
            end
        end
        
        function dNdy = get_dNdy(obj, param_coords)
            % GET_DNDY Computes physical y-derivatives of shape functions (2D only)
            if obj.dim == 1
                error('y-derivative is not defined for a 1D problem.');
            elseif obj.dim == 2
                J = obj.get_jacobian(param_coords);
                
                % --- YOUR CODE HERE FOR 2D ---
                dNdy = zeros(1, size(obj.elem_coords, 1));
            end
        end
    end
end