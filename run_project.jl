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
        "--parallel"
            help = "Use parallel computing"
            action = :store_true
        "--temporal"
            help = "Use time-series train/test split"
            action = :store_true    
        "--rta_only"
            help = "Use RTA trips only"
            action = :store_true    
        "--original"
            help = "Use original trips only"
            action = :store_true
        "--depot"
            help = "Latitude and longitude for depot"
            arg_type = Tuple{Float32, Float32}
            default = (48.440353f0, -123.369201f0)
        "--split"
            help = "Train-test split"
            arg_type = Float32
            default = 0.75
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1
        "--time_limit"
            help = "Time limit for full model in seconds"
            arg_type = Int
            default = 3600
        "--method"
            help = "LP solution algorithm"
            arg_type = Int
            default = -1
        "--delay_cost"
            help = "Cost per passenger-hour of delay"
            arg_type = Float32
        "--orig_thresh"
            help = "Percentage of original trips that must be used"
            arg_type = Float32
        "--stop_time"
            help = "End of the planning horizon"
            arg_type = Int
        "--start_time"
            help = "Start of the planning horizon"
            arg_type = Int
        "--max_vehicles"
            help = "Use the same number of vehicles as the minimum operational cost solution + max_vehicles"
            arg_type = Int
        "--routes"
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
    orig_thresh_str = isnothing(args["orig_thresh"]) ? "" : "_$(args["orig_thresh"])_min_original"
    delay_cost_str = isnothing(args["delay_cost"]) ? "" : "_\$$(args["delay_cost"])"
    max_vehicles_str = isnothing(args["max_vehicles"]) ? "" : "_plus$(args["max_vehicles"])veh"
    data_folder = joinpath(@__DIR__(), "data/objects")
    log_folder = joinpath(@__DIR__(), "outputs")
    filename = "$(routes_str)$(start_time_str)$(stop_time_str)$(
        args["original"] ? "_original" : args["rta_only"] ? "_rta_only" : ""
    )$(
        args["temporal"] ? "_temporal" : ""
    )$(
        isnothing(args["orig_thresh"]) ? "" : orig_thresh_str
    )$(
        isnothing(args["delay_cost"]) ? "" : delay_cost_str
    )$(
        isnothing(args["max_vehicles"]) ? "" : max_vehicles_str
    )"
    log_file = joinpath(log_folder, filename * ".log")
    data_file = joinpath(data_folder, filename * ".jld2")
    if args["debug"]
        println("Parsed args:")
        for (arg, val) in args
            println("  $arg  =>  $val")
        end
        println("filename: $(filename)")

        return 0
    end

    @assert !(args["original"] && args["rta_only"]) "Must select either rta or original trips."

    time_limit = args["time_limit"]

    historical_data = loadHistoricalData("../data/victoria_data.csv")
    trips, _ = loadGTFS("../data/Victoria-GTFS", historical_data)
    
    vic_depot = args["depot"]
    trips_subset = subsetGTFS(trips; routes = args["routes"], start_time = args["start_time"], stop_time = args["stop_time"])
    L_train, L_test = getHistoricalDelays(trips_subset, historical_data; temporal = args["temporal"], split = args["split"], randomSeed = args["seed"])
    r = getHistoricalRidership(trips_subset, historical_data)
    if !isnothing(args["orig_thresh"])
        orig_thresh = round(Int, size(trips_subset, 1) * args["orig_thresh"])
    else
        orig_thresh = nothing
    end
    println("=== Building instance ===")
    @time begin
        if isnothing(args["delay_cost"])
            instance = VSPInstance(trips_subset, r, L_train, L_test; original = args["original"], rta_only = args["rta_only"], depot_loc = vic_depot)
        else
            instance = VSPInstance(trips_subset, r, L_train, L_test; original = args["original"], rta_only = args["rta_only"], depot_loc = vic_depot, delay_cost = args["delay_cost"])
        end
    end
    println("=== Building model ===")
    @time begin
        warmStart = solve!(MCFModel(instance; orig_thresh = orig_thresh,  duplicates = !(args["original"] || args["rta_only"]), timeLimit = 3600))
        if !isnothing(args["max_vehicles"])
            max_vehicles = warmStart.numVehicles + args["max_vehicles"]
        else
            max_vehicles = nothing
        end
        model = VSPModel(instance; max_vehicles = max_vehicles, orig_thresh = orig_thresh, warmStart = warmStart, timeLimit = time_limit, endoftrip = true, logFile = log_file, method = args["method"], silent = false)
        get_current_solution_callback!(model, filename, data_folder)
    end
    println("=== Solving model ===")
    @time begin
        solution = solve!(model; silent = false)
    end
    
    jldsave(
        data_file,
        inst=instance,
        x=solution.x,
        sol_time=solution.solve_time,
        objective_value=solution.objective_value,
        status=termination_status(solution.mod.model)
    )
end

main()