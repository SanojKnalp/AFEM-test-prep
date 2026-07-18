classdef FEMesh
    % FEMESH A simple deal.II-style mesh container for 1D and 2D elements.
    
    properties (SetAccess = private)
        nodes       % [n_nodes x dim] matrix of physical coordinates
        elements    % [n_elements x n_nodes_per_elem] matrix of connectivity (lexicographical)
        dim         % Dimension of the mesh (1 or 2)
    end
    
    methods
        function obj = FEMesh(nodes, elements)
            % CONSTRUCTOR: Initialize the mesh with nodes and elements.
            % Example: mesh = FEMesh(msh.nodes, msh.elements);
            obj.nodes = nodes;
            obj.elements = elements;
            obj.dim = size(nodes, 2);
        end
        
        function node_ids = get_local_nodes(obj, elem_id)
            % GET_LOCAL_NODES Returns the global node IDs for a given element.
            % The output ordering inherently matches the lexicographical 
            % order stored in the 'elements' property.
            %
            % Input:  elem_id  - Integer ID of the element
            % Output: node_ids - Row vector of node IDs
            
            node_ids = obj.elements(elem_id, :);
        end
        
        function coords = get_node_positions(obj, node_ids)
            % GET_NODE_POSITIONS Returns the physical coordinates for given node IDs.
            %
            % Input:  node_ids - Vector of node IDs (e.g., from get_local_nodes)
            % Output: coords   - Matrix of coordinates [length(node_ids) x dim]
            
            coords = obj.nodes(node_ids, :);
        end
        
        function elem_coords = get_element_coordinates(obj, elem_id)
            % GET_ELEMENT_COORDINATES Convenience function to directly fetch 
            % all physical coordinates for a specific element.
            %
            % Input:  elem_id     - Integer ID of the element
            % Output: elem_coords - Matrix of coordinates
            
            node_ids = obj.get_local_nodes(elem_id);
            elem_coords = obj.get_node_positions(node_ids);
        end
    end
end