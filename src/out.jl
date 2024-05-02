using Graphs
using GraphPlot
using JuMP
using Random

"""
    VSPSolution

Solution of `mod`.

# Fields
- `numVehicles::Union{Float64, Int}`: number of vehicles required in the solution.
- `isInt::Bool`: whether the arc decision variables, `x`, are integer or not.
- `x::Union{Matrix{Float64}, Matrix{Int}}`: arc decision variables.
- `s::Vector{Float64}`: propagated delay variables.
- `objective_value::Union{Float64, Vector{Float64}}`: optimal objective value.
- `solve_time::Float64`: computation time for the optimization.
- `mod::Union{VSPModel, FirstStageProblem}`: the optimized model.
"""
struct VSPSolution
    numVehicles::Union{Float64, Int} # number of vehicles used
	isInt::Bool # whether the solution is integer or not
	x::Union{Matrix{Float64}, Matrix{Int}} # link decision values
    s::Vector{Float64} # propagated trip delays
    objective_value::Union{Float64, Vector{Float64}}
    solve_time::Float64
    mod::Union{VSPModel, FirstStageProblem}
end

"""
    solve!(mod::VSPModel)

Optimize the VSP model, `mod`.
"""
function solve!(mod::VSPModel)
    optimize!(mod.model)
    numTrips = mod.inst.n - 1
    x = value.(mod.x)
    s = value.(mod.s)
    isInt = all(isinteger.(x))
    numVehicles = sum(x[1, :])

    @show numTrips
    @show numVehicles
    @show isInt
    @show termination_status(mod.model)
    @show objective_value(mod.model)
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
        $(i): $(round(sol.mod.inst.trips[i, :start_time]; digits = 1)) -  $(round(sol.mod.inst.trips[i, :stop_time]; digits = 1)) \n
        l = $(round(sol.mod.inst.l[i+1]; digits = 1)) \n
        s = $(round(s[i]; digits = 2))
        """
        for i ∈ 1:n-1
    ]

    nodelabel = vcat("depot", triplabels)
    nodesize = vcat(0.5, outdegree(g)[2:end])
    nodelabelsize = vcat(2, [1 for _ ∈ 1:n-1])
    edgelabel = [
        """
        x = $(round(sol.x[src(e), dst(e)]; digits = 1)) \n
        b = $(round(sol.mod.inst.B[src(e), dst(e)]; digits = 1))
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

# function plotVSP_time(sol::Union{VSPSolution, MCFSolution})
#     trips = sol.mod.inst.trips
#     x = convert(Matrix{Bool}, round.(sol.x))
#     schedules = generate_blocks(x)
#     time_plot = plot(;legend=:topright)
#     cmap = range(colorant"red", stop=colorant"blue", length=Int(sol.numVehicles))
#     yflip!(true)

#     for (s, schedule) ∈ enumerate(schedules)
#         this_schedule = schedule .- 1
#         xs = []
#         ys = []
#         dhxs = []
#         dhys = []
#         duration = 0.0
#         for (i, trip) ∈ enumerate(this_schedule)
#             push!(xs, duration, duration+trips[trip, :stop_time]-trips[trip, :start_time], NaN)
#             push!(ys, trips[trip, :start_time], trips[trip, :start_time], NaN)

#             try
#                 last_dur = duration+trips[trip, :stop_time]-trips[trip, :start_time]
#                 duration = trips[this_schedule[i+1], :start_time]-trips[this_schedule[1], :start_time]
#                 push!(dhxs, last_dur, duration, NaN)
#                 push!(dhys, trips[trip, :start_time], trips[this_schedule[i+1], :start_time], NaN)
#             catch
#                 last_dur = duration+trips[trip, :stop_time]-trips[trip, :start_time]
#                 push!(dhxs, last_dur, last_dur)
#                 push!(dhys, trips[trip, :start_time], maximum(trips[:, :stop_time]))
#             end
#         end
#         plot!(
#             xs,
#             ys;
#             lc = cmap[s],
#             lw = 2,
#             label = "vehicle $(s)"
#         )
#         plot!(dhxs, dhys; ls=:dash, lc=:black, la=0.25, label=nothing)
#     end

#     return time_plot
# end

function plotVSP_time(sol::Union{VSPSolution, MCFSolution})
    trips = sol.mod.inst.trips
    x = convert(Matrix{Bool}, round.(sol.x))
    schedules = generate_blocks(x)
    first_trips = [schedule[1] for schedule in schedules]
    schedules = schedules[sortperm(trips[first_trips, :start_time])]
    time_plot = plot(;legend=false)
    vehicle_cmap = range(colorant"red", stop=colorant"blue", length=Int(sol.numVehicles))
    blocks = unique(trips[:, :block_id])
    block_cmap = range(colorant"red", stop=colorant"blue", length=length(blocks))
    yflip!(true)

    counter = 0
    for (s, schedule) ∈ enumerate(schedules)
        this_schedule = schedule
        annot_xs = []
        annot_ys = []
        annots = []
        for (i, trip) ∈ enumerate(this_schedule)
            plot!(
                [trips[trip, :start_time], trips[trip, :stop_time]],
                [counter, counter];
                lc = block_cmap[findfirst(==(trips[trip, :block_id]), blocks)],
                lw = 2
            )
            push!(annot_xs, (trips[trip, :start_time]+trips[trip, :stop_time])/2)
            push!(annot_ys, counter-0.75)
            push!(annots, Plots.text("route $(trips[trip, :route_id])", :black, :center, 4))

            try
                plot!(
                    [trips[trip, :stop_time], trips[this_schedule[i+1], :start_time]],
                    [counter, counter+1];
                    ls = :dash,
                    lc = :black,
                    la = 0.25
                )
            catch
                nothing
            end

            counter += 1
        end
        annotate!(annot_xs, annot_ys, annots)
    end

    return time_plot
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

        push!(schedules, schedule .- 1)
    end

    return schedules
end