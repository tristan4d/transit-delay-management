using Distributed
Distributed.@everywhere begin
    import Pkg
    Pkg.activate("./thesis")
end

using ArgParse

include("./src/utils.jl")
include("./src/data.jl")
include("./src/opt.jl")
include("./src/lagrange.jl")
include("./src/out.jl")
include("./src/metrics.jl")

ArgParse.parse_item(::Type{Vector{Int}}, arg::AbstractString) = map(x -> parse(Int, x), split(strip(arg, ['[', ']']), ","))

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--debug"
            help = "Check arguments"
            action = :store_true
        "--normalize"
            help = "Set mean primary delay to 0"
            action = :store_true
        "--parallel"
            help = "Use parallel computing"
            action = :store_true    
        "--depot"
            help = "Latitude and longitude for depot"
            arg_type = Tuple{Float32, Float32}
            default = (48.440353f0, -123.369201f0)
        "--num_scenarios"
            help = "Number of delay scenarios"
            arg_type = Int
            default = 100
        "--multi"
            help = "Primary delay multiplication factor"
            arg_type = Float32
            default = 1.0
        "--split"
            help = "Train-test split"
            arg_type = Float32
            default = 0.8
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1
        "--time_limit"
            help = "Time limit for full model in seconds"
            arg_type = Int
            default = 3600
        "--stop_time"
            help = "End of the planning horizon"
            arg_type = Int
        "--start_time"
            help = "Start of the planning horizon"
            arg_type = Int
        "routes"
            help = "Routes to be used (e.g., [1,2,3])"
            arg_type = Vector{Int}  # Custom type
            required = true
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    routes_str = join(args["routes"], "-")
    start_time_str = isnothing(args["start_time"]) ? "" : "a$(args["start_time"])"
    stop_time_str = isnothing(args["stop_time"]) ? "" : "b$(args["stop_time"])"
    filename = "victoria_$(routes_str)$(start_time_str)$(stop_time_str)_$(args["num_scenarios"])-$(args["split"])_$(args["multi"])_$(args["seed"]).jld2"
    if args["debug"]
        println("Parsed args:")
        for (arg, val) in args
            println("  $arg  =>  $val")
        end
        println("filename: $(filename)")

        return 0
    end

    full_time = args["time_limit"]
    mean_time = max(60, ceil(Int, full_time / (args["num_scenarios"] * args["split"])))
    easy_time = max(60, ceil(Int, mean_time / 10))

    trips, _ = loadGTFS("../data/Victoria-GTFS")
    historical_data = loadHistoricalData("../data/victoria_ridership_fall_2023.csv"; normalize = args["normalize"])
    
    vic_depot = args["depot"]
    trips_subset = subsetGTFS(trips; routes = args["routes"], start_time = args["start_time"], stop_time = args["stop_time"])
    L_train, L_test = getHistoricalDelays(trips_subset, historical_data, args["num_scenarios"], split = args["split"], multi = args["multi"], randomSeed = args["seed"])
    r = getHistoricalRidership(trips_subset, historical_data)
    instance = VSPInstance(trips_subset, r, L_train, L_test; depot_loc = vic_depot)
    trim_instance = trimInstance(instance; timeLimit=mean_time, parallel=args["parallel"])
    mcf_instance = VSPInstance(trips_subset, r, L_train, L_test; depot_loc = vic_depot, original = true)
    model = MCFModel(mcf_instance; timeLimit = easy_time, duplicates = false)
    println("=== MCF Model ===")
    solution = solve!(model; silent = false)
    mcf_x = solution.x
    mcf_sol_time = solution.solve_time
    mcf_status = termination_status(solution.mod.model)
    model = MCFModel(instance; timeLimit = easy_time)
    println("=== MCFRTA Model ===")
    solution = solve!(model, silent = false)
    mcfrta_x = solution.x
    mcfrta_sol_time = solution.solve_time
    mcfrta_status = termination_status(solution.mod.model)
    model = VSPModel(instance; warmStart = solution, L_train = mean(instance.L_train, dims = 2), timeLimit = mean_time, endoftrip = true)
    println("=== Mean Model ===")
    solution = solve!(model; silent = false)
    mean_x = solution.x
    mean_sol_time = solution.solve_time
    mean_status = termination_status(solution.mod.model)
    model = VSPModel(trim_instance; warmStart = solution, timeLimit = full_time, endoftrip = true)
    println("=== Pruned Model ===")
    solution = solve!(model; silent = false)
    trim_x = solution.x
    trim_sol_time = solution.solve_time
    trim_status = termination_status(solution.mod.model)
    model = VSPModel(instance; warmStart = solution, timeLimit = full_time, endoftrip = true)
    println("=== Full Model ===")
    solution = solve!(model; silent = false)
    del_x = solution.x
    del_sol_time = solution.solve_time
    del_status = termination_status(solution.mod.model)
    println("=== RTA Model ===")
    solution = runTimeAnalysis(trips_subset, r, L_train, L_test; silent = false, timeLimit=easy_time)
    rta_x = solution.x
    rta_sol_time = solution.solve_time
    rta_status = termination_status(solution.mod.model)
    
    folder = joinpath(@__DIR__(), "data/objects")
    jldsave(
        joinpath(folder, filename),
        inst=instance,
        mcf_inst=mcf_instance,
        rta_inst=solution.mod.inst,
        mcf_x=mcf_x,
        mcf_sol_time = mcf_sol_time,
        mcf_status=mcf_status,
        mcfrta_x=mcfrta_x,
        mcfrta_sol_time = mcfrta_sol_time,
        mcfrta_status=mcfrta_status,
        mean_x=mean_x,
        mean_sol_time=mean_sol_time,
        mean_status=mean_status,
        del_x=del_x,
        del_sol_time=del_sol_time,
        del_status=del_status,
        del_x_trim=trim_x,
        del_sol_time_trim=trim_sol_time,
        del_status_trim=trim_status,
        rta_x=rta_x,
        rta_solve_time=rta_sol_time,
        rta_status=rta_status,
        best_cost=getBestPossibleStats(instance; parallel=args["parallel"])
    )
end

main()