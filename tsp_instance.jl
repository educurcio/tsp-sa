module TSPInstance

export Instance
 
mutable struct Instance
    name::String
    comment::String
    type::String
    number_of_nodes::Int64
    edge_weight_type::String
    node_list::Matrix{Float64}
    dist_matrix::Matrix{Int64}
    
    function Instance(filename::String) 
        name, comment, type, number_of_nodes, edge_weight_type, node_list = read_tsp_file(filename)
        dist_matrix = calculate_distance_matrix(edge_weight_type, node_list)
        new(name, comment, type, number_of_nodes, edge_weight_type, node_list, dist_matrix)
    end

    function read_tsp_file(filename::String)
        lines = readlines(filename)

        name = ""
        edge_weight_type = ""
        comment = ""
        type = ""
        dimension = 0
        node_list = zeros(UInt64, 0, 2)

        for line in lines
            # Split the line into key and value based on ":"
            parts = split(line, ":", limit=2)
            key = strip(parts[1])
            value = strip(get(parts, 2, ""))
    
            # Update variables based on the key
            if key == "NAME"
                name = String(value)
            elseif key == "COMMENT"
                comment = String(value)
            elseif key == "TYPE"
                type = String(value)
            elseif key == "DIMENSION"
                dimension = parse(Int, value)
            elseif key == "EDGE_WEIGHT_TYPE"
                edge_weight_type = String(value)
            elseif key == "NODE_COORD_SECTION"
                # Find the starting line for node coordinates
                start_line = findfirst(x -> x == "NODE_COORD_SECTION", lines)
                if start_line !== nothing
                    # Iterate through the coordinates
                    for i in start_line+1:start_line+dimension
                        coord_line = lines[i]
                        node_id, x, y = map(x -> parse(Float64, x), split(coord_line))
                        node_list = vcat(node_list, [x, y]')
                    end
                end
            end
        end
        return name, comment, type, dimension, edge_weight_type, node_list
    end
    
    function calculate_distance_matrix(edge_weight_type::String, node_list::Matrix{Float64})::Matrix{UInt64}
        n_cities = size(node_list, 1)
        dist_matrix = zeros(UInt64, n_cities, n_cities)
        if edge_weight_type == "ATT"
            for i in 1:n_cities
                for j in 1:n_cities
                    dist_matrix[i, j] = ceil((((node_list[i, 1] - node_list[j, 1]) ^ 2 + (node_list[i, 2] - node_list[j, 2]) ^ 2)/10) ^ (1/2))
                end
            end
        else
            for i in 1:n_cities
                for j in 1:n_cities
                    dist_matrix[i, j] = round(((node_list[i, 1] - node_list[j, 1]) ^ 2 + (node_list[i, 2] - node_list[j, 2]) ^ 2) ^ (1/2))
                end
            end
        end
        return dist_matrix
    end
end

end