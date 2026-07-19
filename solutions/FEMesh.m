classdef FEMesh
    % FEMESH A simple deal.II-style mesh container and generator.
    
    properties (SetAccess = private)
        nodes       % [n_nodes x dim] matrix of physical coordinates
        elements    % [n_elements x n_nodes_per_elem] matrix of connectivity
        dim         % Dimension of the mesh (1 or 2)
    end
    
    methods
        function obj = FEMesh(nodes, elements, dim)
            % CONSTRUCTOR: Initialize the mesh manually if you already have the data.
            obj.nodes = nodes;
            obj.elements = elements;
            obj.dim = dim;
        end
        
        function node_ids = get_local_nodes(obj, elem_id)
            % GET_LOCAL_NODES Returns global node IDs for a given element.
            node_ids = obj.elements(elem_id, :);
        end
        
        function coords = get_node_positions(obj, node_ids)
            % GET_NODE_POSITIONS Returns physical coordinates for given node IDs.
            coords = obj.nodes(node_ids, :);
        end
        
        function elem_coords = get_element_coordinates(obj, elem_id)
            % GET_ELEMENT_COORDINATES Fetch all coordinates for a specific element.
            node_ids = obj.get_local_nodes(elem_id);
            elem_coords = obj.get_node_positions(node_ids);
        end

        function [elem_id, xi] = find_element_and_parametric_position(obj, x_target)
            % FIND_ELEMENT_AND_PARAMETRIC_POSITION Finds the element ID and local 
            % parametric coordinate for a given physical position x_target.
            % Note: Currently implemented for 1D meshes.
            
            if obj.dim ~= 1
                error('This method is currently only implemented for 1D meshes.');
            end
            
            n_elements = size(obj.elements, 1);
            elem_id = -1;
            xi = 0;
            
            for e = 1:n_elements
                % Get physical coordinates of the element's nodes
                coords = obj.get_element_coordinates(e);
                x_min = min(coords);
                x_max = max(coords);
                
                % Check if target is within this element (using a tiny tolerance 
                % to catch floating point inaccuracies at element boundaries)
                if x_target >= x_min - 1e-10 && x_target <= x_max + 1e-10
                    elem_id = e;
                    
                    % Map the physical coordinate x to the standard 
                    % parametric domain [-1, 1]
                    xi = 2 * (x_target - x_min) / (x_max - x_min) - 1;
                    
                    % Break early since we found the element
                    return; 
                end
            end
            
            if elem_id == -1
                error('Target point x = %f is outside the mesh domain.', x_target);
            end
        end
    end
    
    methods (Static)
        function obj = generate(dim, n_elem_1d, degree)
            % GENERATE Factory method to create a structured mesh on a unit domain.
            %
            % Inputs:
            %   dim       : 1 or 2
            %   n_elem_1d : Number of elements along one axis
            %   degree    : 1 (linear/bilinear) or 2 (quadratic/biquadratic)
            %
            % Output:
            %   obj       : An instance of FEMesh
            
            if dim == 1
                n_nodes = n_elem_1d * degree + 1;
                nodes = linspace(0, 1, n_nodes)';
                
                elements = zeros(n_elem_1d, degree + 1);
                for i = 1:n_elem_1d
                    start_idx = (i-1)*degree + 1;
                    elements(i, :) = start_idx : (start_idx + degree);
                end
                
            elseif dim == 2
                n_nodes_1d = n_elem_1d * degree + 1;
                
                % Use ndgrid so X changes along rows (lexicographical layout)
                [X, Y] = ndgrid(linspace(0, 1, n_nodes_1d), linspace(0, 1, n_nodes_1d));
                nodes = [X(:), Y(:)];
                
                n_elem = n_elem_1d^2;
                n_loc_nodes = (degree + 1)^2;
                elements = zeros(n_elem, n_loc_nodes);
                
                e = 1;
                for j = 1:n_elem_1d % y-direction
                    for i = 1:n_elem_1d % x-direction
                        % Find the bottom-left node index of the element
                        bl = (j-1)*degree * n_nodes_1d + (i-1)*degree + 1;
                        
                        if degree == 1
                            % Q1: 4 nodes lexicographical
                            elements(e, :) = [bl, bl+1, ...
                                              bl+n_nodes_1d, bl+n_nodes_1d+1];
                        elseif degree == 2
                            % Q2: 9 nodes lexicographical
                            elements(e, :) = [bl, bl+1, bl+2, ...
                                              bl+n_nodes_1d, bl+n_nodes_1d+1, bl+n_nodes_1d+2, ...
                                              bl+2*n_nodes_1d, bl+2*n_nodes_1d+1, bl+2*n_nodes_1d+2];
                        else
                            error('Only degree 1 and 2 are supported for mesh generation.');
                        end
                        e = e + 1;
                    end
                end
            else
                error('Only 1D and 2D mesh generation is supported.');
            end
            
            % Instantiate the object using the standard constructor
            obj = FEMesh(nodes, elements, dim);
        end
    end
end