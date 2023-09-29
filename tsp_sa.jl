module TSPSA

using Random # to use rand
using Base.Threads # For multi-threading
using ..TSPInstance

export initialize_simulated_annealing

@inline function objective_function(solution::Vector{Int64}, dist_matrix::Matrix{Int64})::Int64
    n_cities = size(dist_matrix, 1)
    return sum([dist_matrix[solution[i], solution[i+1]] for i in 1:n_cities])
end
@inline function temperature_schedule(initial_temp::Float64, cooling_factor::Float64)::Float64
    return initial_temp * cooling_factor
end
@inline function acceptance(new_energy::Int64, temperature::Float64)::Bool
    return new_energy <= 0 ? true : @fastmath(exp(Float64(-new_energy) / temperature)) > rand()
end
@inline function calc_new_energy_two_opt(current_solution::Vector{Int64}, dist_matrix::Matrix{Int64}, idx1::Int64, idx2::Int64)::Int64
    return dist_matrix[current_solution[idx1-1],current_solution[idx2]] + dist_matrix[current_solution[idx1],current_solution[idx2+1]] - dist_matrix[current_solution[idx1-1],current_solution[idx1]] - dist_matrix[current_solution[idx2],current_solution[idx2+1]]
end
@inline function sample_two_ints(n_cities::Int64)::Tuple{Int64, Int64}
    idx2 = rand(2:n_cities)
    idx1 = rand(2:n_cities)
    if idx1 > idx2
        return idx2, idx1
    end
    return idx1, idx2
end
function generate_initial_solution(number_of_cities::Int64)::Vector{Int64}
    init_solution = collect(1:number_of_cities)
    Random.shuffle!(init_solution)
    push!(init_solution,init_solution[1])
    
    return init_solution
end
@inline function two_opt(current_solution::Vector{Int64}, dist_matrix::Matrix{Int64}, n_nodes::Int64)::Tuple{Int64, Int64, Int64}
    idx1, idx2 = sample_two_ints(n_nodes)
    return calc_new_energy_two_opt(current_solution, dist_matrix, idx1, idx2), idx1, idx2
end
function simulated_annealing(current_solution::Vector{Int64}, dist_matrix::Matrix{Int64}, temperature::Float64, cooling_factor::Float64, final_temperature::Float64)::Tuple{Vector{Int64}, Int64}
    
    current_energy = objective_function(current_solution, dist_matrix)
    best_solution, best_energy = copy(current_solution), current_energy
    n_nodes = size(dist_matrix, 1)
    while temperature > final_temperature
        energy_delta, idx1, idx2 = two_opt(current_solution, dist_matrix, n_nodes)
        if acceptance(energy_delta,temperature)
            reverse!(current_solution,idx1,idx2)
            current_energy += energy_delta
            if current_energy < best_energy
                best_solution, best_energy = current_solution, current_energy
            end
        end
        temperature = temperature_schedule(temperature,cooling_factor)
    end
    best_energy = objective_function(best_solution, dist_matrix)
    return best_solution, best_energy
end

function initialize_simulated_annealing(instance::TSPInstance.Instance, initial_temp::Float64=10.0, final_temperature::Float64= 0.001, cooling_factor::Float64=0.9999999, number_of_runs::Int64=1, seed::Int64 = 1000)

    Random.seed!(seed);
    number_of_nodes = instance.number_of_nodes
    
    best_energy_vector = zeros(Int64, number_of_runs)
    best_solution_vector = Vector{Vector{Int64}}(undef, number_of_runs)
    for i in 1:number_of_runs
        best_solution_vector[i] = Vector{Int64}(undef,number_of_nodes+1)
    end
    Threads.@threads for i in 1:number_of_runs
        initial_solution = generate_initial_solution(number_of_nodes)
        best_solution_vector[i], best_energy_vector[i] = simulated_annealing(initial_solution, instance.dist_matrix, initial_temp, cooling_factor, final_temperature)
    end
    best_energy_vector_out, best_index = findmin(best_energy_vector)
    return best_energy_vector_out, best_solution_vector[best_index]
end

end