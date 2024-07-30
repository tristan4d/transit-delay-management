using DataFrames
using GeoInterface
using Statistics

struct SolutionStats
    sol::Union{VSPSolution, MCFSolution}
    cost::Float64
    μ::Float64
    μ_10::Float64
    σ::Float64
    σ_10::Float64
    utilization::Float64
    metrics::DataFrame
end

function getSolutionStats(
    sol::Union{VSPSolution, MCFSolution},
    shapes::DataFrame,
    delays::Union{Matrix{Float64}, Nothing} = nothing;
    ridership = nothing
)
    x = convert(Matrix{Bool}, round.(sol.x))
    schedules = generate_blocks(x)
    n = sol.mod.inst.n
    if isnothing(delays)
        L = sol.mod.L_train
        numScenarios = sol.mod.n_train
    else
        L = delays
        numScenarios = size(delays, 2)
    end
    
    B = sol.mod.inst.B
    D = sol.mod.inst.D
    propagated_delays = zeros(Float64, n-1)
    propagated_delay_errs = zeros(Float64, n-1)

    this_s = zeros(Float64, n, numScenarios)
    for scenario in 1:numScenarios
        this_s[:, scenario] = feasibleDelays(sol.x, L[:, scenario], B)
        if !isnothing(ridership)
            this_s[2:end, scenario] .*= ridership
        end
    end
    propagated_delays = vec(mean(this_s, dims=2))[2:end]
    propagated_delay_errs = vec(std(this_s, dims=2))[2:end]
    trips = sol.mod.inst.trips
    metrics = DataFrame(
        [
            Float64[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Vector{Vector{GeoInterface.LineString}}()
        ],
        [
            "duration",
            "num_trips",
            "utilization",
            (isnothing(ridership) ? "propagated_delay" : "propagated_passenger_delay"),
            (isnothing(ridership) ? "propagated_delay_err" : "propagated_passenger_delay_err"),
            "trip_distance",
            "deadhead_distance",
            "geometry"
        ]
    )

    total_duration = 0.0
    total_nis_length = 0.0
    for schedule in schedules
        duration = getBlockLength(schedule, trips)
        total_duration += duration
        num_trips = length(schedule)
        utilization = 1 - notInServiceLength(schedule, trips) / duration
        total_nis_length += notInServiceLength(schedule, trips)
        propagated_delay = mean(propagated_delays[schedule])
        propagated_delay_err = std(propagated_delays[schedule])
        if isnan(propagated_delay_err)
            propagated_delay_err = 0.0
        end
        distance, geometry = getGeometry(schedule, trips, shapes)
        deadhead_distance = getDeadhead(schedule, D)
        push!(metrics, [
            duration,
            num_trips,
            utilization,
            propagated_delay,
            propagated_delay_err,
            distance,
            deadhead_distance,
            geometry
        ])
    end
    
    cost = sum(sol.mod.inst.C .* x) + sum(propagated_delays) * sol.mod.inst.delay_cost + sum(trips[:, :stop_time] .- trips[:, :start_time]) * sol.mod.inst.op_cost
    
    s_10 = ceil(Int, numScenarios/10)
    return SolutionStats(
        sol,
        cost,
        mean(this_s) * 60,
        mean(sort(this_s, dims=2, rev=true)[:, 1:s_10]) * 60,
        std(this_s) * 60,
        std(sort(this_s, dims=2, rev=true)[:, 1:s_10]) * 60,
        1 - total_nis_length / total_duration, 
        metrics
    )
end

function getDeadhead(s::Vector{Int}, D::Matrix{Float64})
    distance = 0.0

    for i in 1:size(s, 1)-1
        distance += D[s[i], s[i+1]]
    end

    return distance
end

function getBlockLength(s::Vector{Int}, trips::DataFrame)
    start = trips[s[1], :start_time]
    stop = trips[s[end], :stop_time]

    return stop - start
end

function notInServiceLength(s::Vector{Int}, trips::DataFrame)
    nis = 0.0

    for i in 1:length(s)-1
        stop = trips[s[i], :stop_time]
        start = trips[s[i+1], :start_time]

        nis += start - stop
    end

    return nis
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
    x::Matrix{Float64},
    l::Vector{Float64},
    B::Matrix{Float64}
)
    n = size(x, 1)
    s = zeros(Float64, n)
    delays = ones(Float64, n)
    iterations = 0

    while any(delays .> eps())
        if iterations > 100
            break
        end

        delays = [x[:, i]' * (s .+ l .- B[:, i]) - s[i] for i ∈ 1:n]

        iterations += 1
        s = max.(0, s .+ delays)
    end

    return s
end

function getGeometry(
    s::Vector{Int},
    trips::DataFrame,
    shapes::DataFrame
)
    coordinates = GeoInterface.LineString[]
    distance = 0.0
    for (t, trip) in enumerate(s)
        shape = shapes[shapes.shape_id .== trips[trip, :shape_id], :shape_pts]

        geom = GeoInterface.LineString(shape[1])
        push!(coordinates, geom)
        try
            dh_start = GeoInterface.Point(trips[trip, :stop_lon], trips[trip, :stop_lat])
            dh_end = GeoInterface.Point(trips[s[t+1], :start_lon], trips[s[t+1], :start_lat])
            push!(coordinates, GeoInterface.LineString([dh_start, dh_end]))
        catch
            nothing
        end

        distance += shapes[shapes.shape_id .== trips[trip, :shape_id], :shape_dist_traveled][1]
    end

    return distance, coordinates
end

function compareSchedules(
    sol_1::Union{VSPSolution, MCFSolution},
    sol_2::Union{VSPSolution, MCFSolution}
)
    x_1 = convert(Matrix{Bool}, round.(sol_1.x))
    x_2 = convert(Matrix{Bool}, round.(sol_2.x))

    return sum(x_1 .& x_2) / sum(x_1 .| x_2)
end