include("./tsp_instance.jl")
include("./tsp_sa.jl") 

using Dates # to use now()
using .TSPInstance
using .TSPSA

path_of_instance = "./instances/dja1436.tsp"
instance = Instance(path_of_instance)

initial_temperature = 10.0
final_temperature = 0.001
cooling_factor = 0.99999999
number_of_runs = 6
seed = 1000

start = now()
energy, best_solution_found = initialize_simulated_annealing(instance, initial_temperature, final_temperature, cooling_factor, number_of_runs, seed)
println(now()-start)
println(energy)