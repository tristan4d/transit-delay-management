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
	x::Union{Matrix{Float64}, Matrix{Int}} # link decision values
    s::Vector{Float64} # propagated trip delays
    s_err::Vector{Float64} # standard deviation of trip delays
    objective_value::Union{Float64, Vector{Float64}}
    solve_time::Float64
    mod::Union{VSPModel, FirstStageProblem}
end

"""
    solve!(mod::VSPModel)

Optimize the VSP model, `mod`.
"""
function solve!(mod::VSPModel; silent=true)
    optimize!(mod.model)
    numTrips = mod.inst.n - 1
    x = value.(mod.x)
    s = value.(mod.s)
    isInt = all(isinteger.(x))
    numVehicles = sum(x[1, :])

    if !silent
        @show numTrips
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

"""
    solve!(mod::FirstStageProblem)

Optimize the first-stage problem model, `mod`.
"""
function solve!(mod::FirstStageProblem)
    optimize!(mod.model)
    numTrips = mod.inst.n - 1
    x = value.(mod.x)
    s = feasibleDelays(x, mod.inst.l, mod.inst.B)
    isInt = all(isinteger.(x))
    numVehicles = sum(x[1, :])

    @show numTrips
    @show numVehicles
    @show isInt
    @show termination_status(mod.model)
    @show objective_value(mod.model)
    @show value(mod.q)
    @show solve_time(mod.model)

    return VSPSolution(
        numVehicles,
        all(isinteger.(x)),
        x,
        s,
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
	x::Union{Matrix{Float64}, Matrix{Int}} # link decision values
    objective_value::Float64
    solve_time::Float64
    mod::MCFModel
end

"""
    solve!(mod::MCFModel)

Optimize the min-cost flow model, `mod`.
"""
function solve!(mod::MCFModel)
    optimize!(mod.model)
    x = value.(mod.x)
    numVehicles = sum(x[1, :])

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

Colors indicate the block in the original GTFS files.  The number indicates the route
number.
"""
function plotVSP_time(
    sol::Union{VSPSolution, MCFSolution};
    delays = nothing,
    ridership = nothing
    )
    trips = sol.mod.inst.trips
    B = sol.mod.inst.B
    D = sol.mod.inst.D
    x = convert(Matrix{Bool}, round.(sol.x))
    s = nothing

    if isnothing(delays)
        s = sol.s
    else
        this_s = zeros(Float64, size(delays))
        for scenario in 1:size(delays, 2)
            this_s[:, scenario] = feasibleDelays(sol.x, delays[:, scenario], B)
        end
        s = vec(mean(this_s, dims=2))[2:end]
    end
    if !isnothing(ridership)
        s = s .* ridership
    end

    schedules = generate_blocks(x)
    delay_cmap = reverse(ColorSchemes.roma)
    time_plot = plot(;
        xlabel="time of day (hrs)",
        ylabel="vehicle schedule",
        yticks=0:sol.numVehicles+1,
        legend=false,
        colorbar=true
    )
    # blocks = unique(trips[:, :block_id])
    # block_cmap = range(colorant"yellow", stop=colorant"blue", length=length(blocks))
    yflip!(true)

    counter = 1
    for schedule ∈ schedules
        this_schedule = schedule
        annot_xs = []
        annot_ys = []
        annots = []

        plot!(
            [
                trips[this_schedule[1], :start_time]-D[1, this_schedule[1]],
                trips[this_schedule[1], :start_time]
            ],
            [counter, counter];
            ls = :solid,
            lc = :black,
            la = 0.75,
            lw = 2
        )
        plot!(
            [
                trips[this_schedule[end], :stop_time],
                trips[this_schedule[end], :stop_time]+D[this_schedule[end], 1]
            ],
            [counter, counter];
            ls = :solid,
            lc = :black,
            la = 0.75,
            lw = 2
        )
        for (i, trip) ∈ enumerate(this_schedule)
            plot!(
                [trips[trip, :start_time], trips[trip, :stop_time]],
                [counter, counter];
                # lc = block_cmap[findfirst(==(trips[trip, :block_id]), blocks)],
                lc = get(delay_cmap, (s[trip]-minimum(s))/(maximum(s)-minimum(s))),
                lw = 10
            )
            push!(annot_xs, (trips[trip, :start_time]+trips[trip, :stop_time])/2)
            push!(annot_ys, counter-0.25)
            push!(annots, Plots.text(trips[trip, :route_id], :black, :center, 4))

            try
                start = trips[trip, :stop_time]
                stop = trips[this_schedule[i+1], :start_time]
                if stop - start < 3
                    deadhead = D[trip, this_schedule[i+1]]
                    plot!(
                        [start, start+deadhead],
                        [counter, counter];
                        ls = :solid,
                        lc = :black,
                        la = 0.75,
                        lw = 2
                    )
                    plot!(
                        [start+deadhead, stop],
                        [counter, counter];
                        ls = :solid,
                        lc = :black,
                        la = 0.75
                    )
                else
                    plot!(
                        [start, stop],
                        [counter, counter];
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

    ylims!(0, counter)

    gr()
    l = @layout [a{0.95w} b]
    cmap = cgrad(delay_cmap)
    p2 = heatmap(
        rand(2,2),
        clims=(minimum(s)*60, maximum(s)*60),
        framestyle=:none,
        c=cmap,
        cbar=true,
        lims=(-1,0),
        colorbar_title=(isnothing(ridership) ? "delay (mins)" : "passenger delay (passenger ⋅ mins)")
    )

    return plot(time_plot, p2, layout=l)
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

function runTimeAnalysis(sol::MCFSolution, delays::Matrix{Float64}; percentile = nothing)
    trips = copy(sol.mod.inst.trips)
    if isnothing(percentile)
        delay_cost = sol.mod.inst.delay_cost
        op_cost = sol.mod.inst.op_cost
        r = sol.mod.inst.r
        percentile = max.(1 .- op_cost ./ delay_cost ./ r, 0)
        pushfirst!(percentile, 0)
    end

    q = quantile.(eachrow(delays), percentile)
    trips.stop_time .+= q[2:end]

    instance = VSPInstance(trips)

    rta_model = MCFModel(instance)
    rta_solution = solve!(rta_model)

    return rta_solution, q
end