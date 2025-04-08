using Distributed
using Graphs
using GraphPlot
using Plots
using Leaflet
using GeoInterface
using Colors
using ColorSchemes
using JuMP
using Random
using Statistics
using StatsBase

Distributed.@everywhere begin
    """
        VSPSolution

    Solution of `mod`.

    # Fields
    - `numVehicles::Union{Float64, Int}`: number of vehicles required in the solution.
    - `isInt::Bool`: whether the arc decision variables, `x`, are integer or not.
    - `x::Union{Matrix{Float64}, Matrix{Int}}`: arc decision variables.
    - `s::Vector{Float64}`: mean trip delays across all scenarios.
    - `s_err::Vector{Float64}`: standard deviation of trip delays across all scenarios.
    - `objective_value::Union{Float64, Vector{Float64}}`: optimal objective value.
    - `solve_time::Float64`: computation time for the optimization.
    - `mod::Union{VSPModel, FirstStageProblem}`: the optimized model.
    """
    struct VSPSolution
        numVehicles::Union{Float64, Int} # number of vehicles used
        isInt::Bool # whether the solution is integer or not
        x::Matrix{Bool} # link decision values
        s::Vector{Float64} # propagated trip delays
        s_err::Vector{Float64} # standard deviation of trip delays
        objective_value::Union{Float64, Vector{Float64}}
        solve_time::Float64
        mod::Union{VSPModel, FirstStageProblem}
    end
end

Distributed.@everywhere begin
    """
        solve!(mod::VSPModel)

    Optimize the VSP model, `mod`.
    """
    function solve!(mod::VSPModel; silent=true)
        optimize!(mod.model)
        numTrips = mod.inst.n - 1
        numScenarios = size(mod.inst.L_train, 2)
        x = zeros(Bool, mod.inst.n, mod.inst.n)
        for (i, j) in Tuple.(findall(mod.inst.G))
            x[i, j] = round(value(mod.x[i, j]))
        end
        s = value.(mod.s)
        isInt = all(isinteger.(x))
        numVehicles = sum(x[1, :])

        if !silent
            @show numTrips
            @show numScenarios
            @show numVehicles
            @show isInt
            @show termination_status(mod.model)
            @show objective_value(mod.model)
            @show solve_time(mod.model)
        end

        return VSPSolution(
            numVehicles,
            all(isinteger.(x)),
            x,
            vec(mean(s, dims=2)),
            vec(std(s, dims=2)),
            objective_value(mod.model),
            solve_time(mod.model),
            mod
        )
    end
end

"""
    solve!(mod::FirstStageProblem)

Optimize the first-stage problem model, `mod`.
"""
function solve!(mod::FirstStageProblem)
    optimize!(mod.model)
    numTrips = mod.inst.n - 1
    L_train = mod.inst.L_train
    B = mod.inst.B
    r = mod.inst.r
    x = value.(mod.x)
    s = zeros(Float64, size(L_train))
    for i in 1:size(L_train, 2)
        s[:, i] .= feasibleDelays(x, L_train[:, i], B, endoftrip=true)
        s[2:end, i] .*= r
    end
    isInt = all(isinteger.(x))
    numVehicles = sum(x[1, :])

    @show numTrips
    @show numVehicles
    @show isInt
    @show termination_status(mod.model)
    @show objective_value(mod.model)
    # @show value(mod.q)
    @show solve_time(mod.model)

    return VSPSolution(
        numVehicles,
        all(isinteger.(x)),
        x,
        vec(mean(s, dims=2)),
        vec(std(s, dims=2)),
        objective_value(mod.model),
        solve_time(mod.model),
        mod
    )
end

"""
    MCFSolution

Solution of `mod`.

# Fields
- `numVehicles::Union{Float64, Int}`: number of vehicles required in the solution.
- `x::Union{Matrix{Float64}, Matrix{Int}}`: arc decision variables.
- `objective_value::Union{Float64, Vector{Float64}}`: optimal objective value.
- `solve_time::Float64`: computation time for the optimization.
- `mod::Union{VSPModel, FirstStageProblem}`: the optimized model.
"""
struct MCFSolution
    numVehicles::Union{Float64, Int} # number of vehicles used
	x::Matrix{Bool} # link decision values
    objective_value::Float64
    solve_time::Float64
    mod::MCFModel
end

"""
    solve!(mod::MCFModel)

Optimize the min-cost flow model, `mod`.
"""
function solve!(mod::MCFModel; silent = true)
    optimize!(mod.model)
    numTrips = mod.inst.n - 1
    x = zeros(Bool, mod.inst.n, mod.inst.n)
    for (i, j) in Tuple.(findall(mod.inst.G))
        x[i, j] = round(value(mod.x[i, j]))
    end
    numVehicles = sum(x[1, :])

    if !silent
        @show numTrips
        @show numVehicles
        @show termination_status(mod.model)
        @show objective_value(mod.model)
        @show solve_time(mod.model)
    end

    return MCFSolution(
        numVehicles,
        x,
        objective_value(mod.model),
        solve_time(mod.model),
        mod
    )
end

# (NOT USING)
function solve!(mod::DelayModel)
    optimize!(mod.model)
    s = value.(mod.s)

    return s
end

"""
    plotVSP(sol::Union{VSPSolution, MCFSolution})

Plot `sol` as a graph.
"""
function plotVSP(sol::Union{VSPSolution, MCFSolution})
    g = SimpleDiGraph(sol.x)
    n = nv(g)
    s = zeros(Float64, n)

    try
        s = sol.s
    catch
        nothing
    end

    triplabels = [
        """
        $(i): $(round(sol.mod.inst.trips[i, :start_time]; digits=1)) -  $(round(sol.mod.inst.trips[i, :stop_time]; digits=1)) \n
        l = $(round(mean(sol.mod.inst.l[i]); digits=1)) +/- $(round(std(sol.mod.inst.l[i]); digits=1))\n
        s = $(round(s[i]; digits = 2))
        """
        for i ∈ 1:n-1
    ]

    nodelabel = vcat("depot", triplabels)
    nodesize = vcat(0.5, outdegree(g)[2:end])
    nodelabelsize = vcat(2, [1 for _ ∈ 1:n-1])
    edgelabel = [
        """
        x = $(round(sol.x[src(e), dst(e)]; digits=1)) \n
        b = $(round(sol.mod.inst.B[src(e), dst(e)]; digits=1))
        """
        for (i, e) ∈ enumerate(edges(g))
    ]
    edgelabelsize = [1 + sol.x[src(e), dst(e)] for (i, e) ∈ enumerate(edges(g))]
    
    return gplot(
        g,
        nodelabelsize = nodelabelsize,
        nodesize = nodesize,
        nodelabel = nodelabel,
        nodelabeldist = 1,
        nodelabelangleoffset = 3*pi/2,
        edgelabelsize = edgelabelsize,
        edgelabeldistx = 0,
        edgelabeldisty = 0,
        edgelabel = edgelabel
    )
end

"""
    plotVSP_time(
    sol::Union{VSPSolution, MCFSolution}[,
    delays = nothing,
    ridership = nothing
    ]
    )

Plot the vehicle schedules in `sol` over time.

Colors indicate the secondary passenger delay cost.  The number indicates the route
number.
"""
function plotVSP_time(
    x::Union{Matrix{Float64}, Matrix{Bool}},
    instance::VSPInstance;
    plot_size = (800, 600),
    endoftrip = true,
    show_rta = false,
    delays = nothing,
    ridership = nothing,
    clims = nothing,
    primary = false,
    sim_no_rta = false,
    y_skip = 1,
    left_margin = :match,
    show_dir = false,
    show_trip = false,
    show_cbar = true,
    show_ylabel = true,
    show_xlabel = true,
    title = ""
    )
    trips = copy(instance.trips)
    rta_mask = instance.rta_mask
    n_original =  size(trips, 1) - sum(rta_mask)
    B = instance.B
    D = instance.D
    x = convert(Matrix{Bool}, round.(x))

    if isnothing(delays)
        delays = instance.L_test
    end
    if isnothing(ridership)
        ridership = instance.r
    end

    if primary
        s = vec(mean(delays, dims=2))[2:end]
        if sim_no_rta
            s[n_original+1:end] .= s[rta_mask]
            trips[n_original+1:end, :] .= trips[rta_mask, :]
        end
    else
        this_s = zeros(Float64, size(delays))
        for scenario in 1:size(delays, 2)
            this_s[:, scenario] = feasibleDelays(x, delays[:, scenario], B, endoftrip=endoftrip)
        end
        s = vec(mean(this_s, dims=2))[2:end]
    end
    s = s .* ridership

    schedules = generate_blocks(x)
    # delay_cmap = reverse(ColorSchemes.roma)
    delay_cmap = ColorSchemes.Reds

    if isnothing(clims)
        clims = (0, max(maximum(s), 1))
    end
    time_plot = plot(;
        xlabel=show_xlabel ? "time of day (hours)" : "",
        ylabel=show_ylabel ? "vehicle schedule" : "",
        yticks=(1:length(schedules), [i % y_skip == 0 ? i : "" for i in 1:length(schedules)]),
        ylims=(0, length(schedules)+1),
        left_margin=left_margin,
        legend=false,
        # legend=:outertopright,
        colorbar=true,
        size=plot_size,
        title=title
    )
    yflip!(true)

    rectangle(w, h, x, y) = Plots.Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    counter = 1
    num_rta = 0
    for schedule ∈ schedules
        this_schedule = schedule
        annot_xs = []
        annot_ys = []
        annots = []

        plot!(
            [
                trips[this_schedule[1], :start_time]-D[1, this_schedule[1]+1],
                trips[this_schedule[1], :start_time]
            ],
            [counter, counter];
            label="",
            ls = :solid,
            lc = :black,
            la = 0.75,
            lw = 2
        )
        plot!(
            [
                trips[this_schedule[end], :stop_time],
                trips[this_schedule[end], :stop_time]+D[this_schedule[end]+1, 1]
            ],
            [counter, counter];
            label="",
            ls = :solid,
            lc = :black,
            la = 0.75,
            lw = 2
        )
        for (i, trip) ∈ enumerate(this_schedule)
            if trip > n_original
                rta = true
                num_rta += 1
                paired_trip = findall(x -> x, rta_mask)[trip-n_original]
            else
                rta = false
            end
            if show_rta
                if rta
                    plot!(
                        rectangle(trips[trip, :stop_time]-trips[trip, :start_time], 0.6, trips[trip, :start_time], counter-0.3),
                        label="",
                        lc=:blue,
                        c=:grey
                    )
                    plot!(
                        rectangle(trips[paired_trip, :stop_time]-trips[trip, :start_time], 0.6, trips[trip, :start_time], counter-0.3),
                        label="",
                        lc=:blue,
                        c=get(delay_cmap, (s[trip]-clims[1])/(clims[2]-clims[1]))
                    )
                    push!(annot_xs, (trips[trip, :start_time]+trips[paired_trip, :stop_time])/2)
                else
                    plot!(
                        rectangle(trips[trip, :stop_time]-trips[trip, :start_time], 0.6, trips[trip, :start_time], counter-0.3),
                        label="",
                        lc=:black,
                        c=get(delay_cmap, (s[trip]-clims[1])/(clims[2]-clims[1]))
                    )
                    push!(annot_xs, (trips[trip, :start_time]+trips[trip, :stop_time])/2)
                end
            else
                plot!(
                    rectangle(trips[trip, :stop_time]-trips[trip, :start_time], 0.6, trips[trip, :start_time], counter-0.3),
                    label="",
                    lc=:black,
                    c=get(delay_cmap, (s[trip]-clims[1])/(clims[2]-clims[1]))
                )
                push!(annot_xs, (trips[trip, :start_time]+trips[trip, :stop_time])/2)
            end
            push!(annot_ys, counter)
            push!(
                annots,
                Plots.text(
                    "$(show_trip ? rta ? paired_trip : trip : trips[trip, :route_id])$(show_dir ? "\n" * trips[trip, :direction_code] : "")",
                    6,
                    :hcenter,
                    :vcenter,
                    (s[trip]-clims[1])/(clims[2]-clims[1]) < 0.5 ? :black : :white
                    )
            )

            try
                start = trips[trip, :stop_time]
                stop = trips[this_schedule[i+1], :start_time]
                if stop - start < 3
                    deadhead = D[trip+1, this_schedule[i+1]+1]
                    plot!(
                        [start, start+deadhead],
                        [counter, counter];
                        label="",
                        ls = :solid,
                        lc = :black,
                        la = 0.75,
                        lw = 2
                    )
                    plot!(
                        [start+deadhead, stop],
                        [counter, counter];
                        label="",
                        ls = :solid,
                        lc = :black,
                        la = 0.75
                    )
                else
                    plot!(
                        [start, stop],
                        [counter, counter];
                        label="",
                        ls = :dash,
                        lc = :black,
                        la = 0.75
                    )
                end
            catch
                nothing
            end
        end

        counter += 1
        annotate!(annot_xs, annot_ys, annots)
    end

    # ylims!(time_plot, (0, counter))
    # yticks!(time_plot, 1:counter-1)

    gr()
    cmap = cgrad(delay_cmap)
    l = @layout [
        a{0.95w} b{0.05w}
    ]
    p2 = heatmap(
        rand(2,2),
        clims=(clims[1], clims[2]),
        framestyle=:none,
        c=cmap,
        cbar=true,
        lims=(-1,0),
        colorbar_title=(primary ? (sim_no_rta ? "pre-RTA primary delay (passenger-mins)" : "primary delay (passenger-mins)") : "passenger delay (passenger-hours)")
    )
    show_rta && println(num_rta / n_original)
    return show_cbar ? plot(time_plot, p2, layout=l) : time_plot
end

"""
    plotVSP(inst::VSPInstance)

Plot `inst` as a graph.
"""
function plotVSP(inst::VSPInstance)
    g = SimpleDiGraph(inst.G)
    n = nv(g)

    triplabels = [
        """
        $(i): $(round(inst.trips[i, :start_time]; digits = 1)) -  $(round(inst.trips[i, :stop_time]; digits = 1)) \n
        l = $(round(inst.l[i+1]; digits = 1)) \n
        s = ?
        """
        for i ∈ 1:n-1
    ]

    nodelabel = vcat("depot", triplabels)
    nodelabelsize = vcat(2, [1 for _ ∈ 1:n-1])
    
    return gplot(
        g,
        nodelabelsize = nodelabelsize,
        nodelabel = nodelabel,
        nodelabeldist = 1,
        nodelabelangleoffset = 3*pi/2,
    )
end

function plotVSP_map(metrics::DataFrame; schedule = nothing)
    lines = Leaflet.Layer[]
    endpoints = Leaflet.Layer[]
    # n = size(metrics, 1)
    # block_cmap = range(colorant"yellow", stop=colorant"blue", length=n)
    block_cmap = palette(:roma50)
    lats = []
    lons = []
    for (i, geom) in enumerate(metrics.geometry)
        if !isnothing(schedule) && i != schedule
            continue
        end
        color = "#" * hex(block_cmap[i])
        translation_vector = rand(0:1, 2) / 5000
        for (j, linestring) in enumerate(geom)
            coords = GeoInterface.coordinates(linestring)
            append!(lons, [coord[1] for coord in coords])
            append!(lats, [coord[2] for coord in coords])
            translated_linestring = GeoInterface.LineString(map(x -> x + translation_vector, coords))
            push!(lines, Leaflet.Layer(
                translated_linestring;
                color=j%2==0 ? "gray" : color,
                opacity=j%2==0 ? 0.5 : 1,
                border_width=2
            ))
            j%2==1 && push!(endpoints, Leaflet.Layer(
                coords[1]+translation_vector;
                color=:white,
                opacity=1,
                marker_size=4
            ))
            j%2==1 && push!(endpoints, Leaflet.Layer(
                coords[end]+translation_vector;
                color=:black,
                opacity=1,
                marker_size=4
            ))
        end
    end

    provider = Leaflet.CARTO()
    m = Leaflet.Map(;
        layers=vcat(lines, endpoints),
        provider=provider,
        zoom=12,
        center=[mean(lats), mean(lons)]
    )

    return m
end

"""
    generate_blocks(x::Matrix{Bool})

Determine the vehicle blocks (schedules) for each vehicle in the solution associated
with `x`.
"""
function generate_blocks(x::Matrix{Bool})
    pull_trips = findall(x[1, :])
    schedules = []

    for trip in pull_trips
        schedule = [trip]
        next_trip = trip

        while findfirst(x[next_trip, :]) > 1
            next_trip = findfirst(x[next_trip, :])
            push!(schedule, next_trip)
        end

        # '.-1' to correct for the depot
        push!(schedules, schedule .- 1)
    end

    return schedules
end

"""
    feasibleDelays(
        x::Matrix{Float64},
        l::Vector{Float64},
        B::Matrix{Float64}
    )

Calculate the propagated delay at each trip given arc decisions `x`, expected delays
`l`, and buffer times `B`.
"""
function feasibleDelays(
    x::Union{Matrix{Float64}, Matrix{Bool}},
    l::Vector{Float64},
    B::Matrix{Float64};
    endoftrip = true
)
    n = size(x, 1)
    s = zeros(Float64, n)
    delays = ones(Float64, n)
    delays[1] = 0
    iterations = 0

    while any(delays .> eps())
        if iterations > 100
            println("max iterations exceeded")
            break
        end

        # delays = [x[:, i+1]' * (s .+ l .- B[:, i]) - s[i] for i ∈ 1:n]
        for i in 1:n-1
            delays[i+1] = x[:, i+1]' * (s .+ l .- B[:, i+1]) - s[i+1]
        end

        iterations += 1
        s = max.(0, s .+ delays)
    end

    if endoftrip
        s = max.(0, s .+ l) .* (sum(x, dims=2) .> 0)
    end

    return s
end

function runTimeAnalysis(
    trips::DataFrame,
    r::Vector{Float64},
    L_train::Matrix{Float64},
    L_test::Matrix{Float64};
    silent = true,
    timeLimit = 600
)
    instance = VSPInstance(copy(trips), r, copy(L_train), copy(L_test), rta_only = true)

    rta_model = MCFModel(instance, duplicates = false, timeLimit=timeLimit)
    rta_solution = solve!(rta_model; silent = silent)

    return rta_solution
end

function trimInstance(
    inst::VSPInstance;
    num_scenarios = nothing,
    parallel = false,
    timeLimit = 3600
)
    G = zeros(Bool, size(inst.G))
    sol = nothing
    m = size(inst.L_train, 2)
    if isnothing(num_scenarios)
        col_means = mean.(eachcol(inst.L_train))
        sorted_indices = sortperm(col_means)
        scenarios = [inst.L_train[:, sorted_indices[1]], inst.L_train[:, sorted_indices[1]]]
        push!(scenarios, mean.(eachrow(inst.L_train)))
    elseif num_scenarios <= m
        indices = sample(1:m, num_scenarios, replace = false)
        scenarios = [inst.L_train[:, i] for i in indices]
    else
        indices = sample(1:m, num_scenarios, replace = true)
        scenarios = [inst.L_train[:, i] for i in indices]
    end
    if parallel
        sols = Distributed.pmap(s -> solve!(VSPModel(inst; L_train=s, timeLimit = timeLimit)), scenarios)

        for sol in sols
            this_x = convert(Matrix{Bool}, round.(sol.x))
            G .|= this_x
        end
    else
        for s in scenarios
            if isnothing(sol)
                mod = VSPModel(inst; L_train=s, timeLimit = timeLimit)
            else
                mod = VSPModel(inst; warmStart=sol, L_train=s, timeLimit = timeLimit)
            end
            sol = solve!(mod)
            this_x = convert(Matrix{Bool}, round.(sol.x))
            G = G .| this_x
            GC.gc()
        end
    end

    this_inst = deepcopy(inst)
    this_inst.G .*= G

    return this_inst
end

function getBestPossibleStats(
    inst::VSPInstance;
    delays = nothing,
    parallel = false,
    timeLimit = 3600,
    warm_start = nothing,
    max_vehicles = nothing
)
    if parallel
        args = eachcol(train ? inst.L_train : inst.L_test)
        obj = Distributed.pmap(arg -> solve!(VSPModel(inst; L_train=arg, timeLimit=timeLimit, max_vehicles=max_vehicles), silent=true).objective_value, args)
    else
        sol = nothing
        obj = []
        for l in eachcol(isnothing(delays) ? inst.L_test : delays)
            if isnothing(warm_start)
                mod = VSPModel(inst; warmStart=sol, L_train=l, timeLimit=timeLimit, endoftrip=true, method=1, max_vehicles=max_vehicles)
            else
                mod = VSPModel(inst; warmStart=warm_start, L_train=l, timeLimit=timeLimit, endoftrip=true, method=1, max_vehicles=max_vehicles)
            end
            sol = solve!(mod; silent=true)
            stats = getSolutionStats(sol.x, inst; delays=l)
            push!(obj, stats.cost)
            GC.gc()
        end
    end
        
    return mean(obj), std(obj), quantile(obj, 0.25), quantile(obj, 0.75)
end

function run_and_save(
    trips,
    routes,
    historical_data;
    temporal = false,
    random_seed = 1,
    time_limit = 7200,
    split = 0.5,
    new_mean = nothing,
    new_std = nothing,
    test_mean = nothing,
    test_std = nothing,
    overlap = 0.0,
    depot = (48.440353, -123.369201),
    max_dist = "mean",
    max_layover = 1.0,
    depot_return_wait = 1.0,
    delay_cost = 36.54,
    veh_cost = 806.10,
    add_vehicles = nothing,
    silent = true,
    method = 1,
    do_rtv = false
)
    trips_subset = subsetGTFS(trips; routes=routes)
    L_train, L_test = getHistoricalDelays(trips_subset, historical_data; temporal=temporal, randomSeed=random_seed, split=split, new_mean=new_mean, new_std=new_std, new_test_mean=test_mean, new_test_std=test_std, overlap=overlap)
    r = getHistoricalRidership(trips_subset, historical_data)
    instance = VSPInstance(trips_subset, r, L_train, L_test; depot_loc=depot, max_dist=max_dist, max_layover=max_layover, depot_return_wait=depot_return_wait, delay_cost=delay_cost, veh_cost=veh_cost)
    if !isnothing(add_vehicles)
        mcf_model = MCFModel(instance; duplicates=true)
        mcf_solution = solve!(mcf_model; silent=silent)
        max_vehicles = mcf_solution.numVehicles + add_vehicles
    else
        max_vehicles = nothing
    end
    mean_model = VSPModel(instance; L_train=mean(instance.L_train, dims=2), timeLimit=time_limit, endoftrip=true, method=method, max_vehicles=max_vehicles) 
    mean_solution = solve!(mean_model; silent=silent)
    full_model = VSPModel(instance; warmStart=mean_solution, timeLimit=time_limit, endoftrip=true, method=method, max_vehicles=max_vehicles)
    full_solution = solve!(full_model; silent=silent)
    if do_rtv
        rtv_instance =  VSPInstance(trips_subset, r, L_train, L_test; depot_loc=depot, max_dist=max_dist, max_layover=max_layover, depot_return_wait=depot_return_wait, delay_cost=delay_cost, rta_only=true)
        rtv_model = VSPModel(rtv_instance; timeLimit=time_limit, endoftrip=true, method=method)
        rtv_solution = solve!(rtv_model; silent=silent)
    end
    best_cost, best_cost_std, best_cost_lo, best_cost_hi = getBestPossibleStats(instance, warm_start=full_solution, max_vehicles=max_vehicles);

    folder = joinpath(@__DIR__(), "..\\data\\objects")
    routes_str = join(routes, "-")
    filename = "$(routes_str)_veh_$(veh_cost).jld2"
    jldsave(
        joinpath(folder, filename),
        inst=instance,
        rtv_inst=do_rtv ? rtv_instance : nothing,
        mean_x=mean_solution.x,
        mean_sol_time=mean_solution.solve_time,
        mean_status=termination_status(mean_solution.mod.model),
        mean_obj=mean_solution.objective_value,
        full_x=full_solution.x,
        full_sol_time=full_solution.solve_time,
        full_status=termination_status(full_solution.mod.model),
        full_obj=full_solution.objective_value,
        rtv_x=do_rtv ? rtv_solution.x : nothing,
        rtv_sol_time=do_rtv ? rtv_solution.solve_time : nothing,
        rtv_status=do_rtv ? termination_status(rtv_solution.mod.model) : nothing,
        rtv_obj=do_rtv ? rtv_solution.objective_value : nothing,
        best_cost=best_cost,
        best_cost_std=best_cost_std,
        best_cost_lo=best_cost_lo,
        best_cost_hi=best_cost_hi
    )
end