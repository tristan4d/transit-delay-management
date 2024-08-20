using DataFrames
using GeoInterface
using Statistics

"""
    SolutionStats

Solution statistics and metrics for minimum cost flow or delay-aware models.

# Fields
- `sol::Union{VSPSolution, MCFSolution}`: a solution of a VSP instance.
- `cost::Float64`: cost, in monetary units, of the solution.
- `cost_err::Float64`: cost confidence interval, in monetary units, of the solution.
- `vehicle_cost::Float64`: cost, in monetary units, of the fleet.
- `link_cost::Float64`: cost, in monetary units, of the chosen links.
- `passenger_cost::Float64`: cost, in monetary units, of the total passenger delay.
- `passenger_cost_err::Float64`: cost confidence interval, in monetary units, of the total passenger delay.
- `service_cost::Float64`: cost, in monetary units, of the service delivered.
- `service_cost_err::Float64`: cost confidence interval, in monetary units, of the service delivered.
- `μ::Float64`: the mean (passenger) delay per trip.
- `μ_10::Float64`: the mean (passenger) delay per trip of the worst 10% of instances.
- `σ::Float64`: the standard deviation of (passenger) delay per trip.
- `σ_10::Float64`: the standard deviation of (passenger) delay per trip of the worst 10% of instances.
- `utilization::Float64`: the percent of time spent moving passengers.
- `deadhead::Float64`: the distance, in minutes, of deadheading in the solution.
- `metrics::DataFrame`: trip-level metrics.
"""
struct SolutionStats
    sol::Union{VSPSolution, MCFSolution}
    cost::Float64
    cost_err::Float64
    vehicle_cost::Float64
    link_cost::Float64
    passenger_cost::Float64
    passenger_cost_err::Float64
    service_cost::Float64
    service_cost_err::Float64
    μ::Float64
    μ_10::Float64
    σ::Float64
    σ_10::Float64
    utilization::Float64
    deadhead::Float64
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
    op_cost = sol.mod.inst.op_cost
    delay_cost = sol.mod.inst.delay_cost
    veh_cost = sol.mod.inst.veh_cost
    B = sol.mod.inst.B
    C = sol.mod.inst.C
    D = sol.mod.inst.D
    V = sol.mod.inst.V
    trips = sol.mod.inst.trips
    propagated_delays = zeros(Float64, n-1)
    propagated_delay_errs = zeros(Float64, n-1)

    if isnothing(delays)
        L = sol.mod.L_train
        numScenarios = sol.mod.n_train
    else
        L = delays
        numScenarios = size(delays, 2)
    end
    this_s = zeros(Float64, n, numScenarios)

    for scenario in 1:numScenarios
        this_s[:, scenario] = feasibleDelays(sol.x, L[:, scenario], B)
        if !isnothing(ridership)
            this_s[2:end, scenario] .*= ridership
        end
    end

    propagated_delays = vec(mean(this_s, dims=2))[2:end]
    propagated_delay_errs = vec(std(this_s, dims=2))[2:end]
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
            "deadhead",
            "geometry"
        ]
    )

    total_duration = 0.0
    total_nis_length = 0.0
    total_deadhead = 0.0
    for schedule in schedules
        duration = getBlockLength(schedule, trips, D)
        total_duration += duration
        num_trips = length(schedule)
        utilization = 1 - notInServiceLength(schedule, trips, D) / duration
        total_nis_length += notInServiceLength(schedule, trips, D)
        propagated_delay = mean(propagated_delays[schedule])
        propagated_delay_err = std(propagated_delays[schedule])
        if isnan(propagated_delay_err)
            propagated_delay_err = 0.0
        end
        distance, geometry = getGeometry(schedule, trips, shapes)
        deadhead = getDeadhead(schedule, D)
        total_deadhead += deadhead
        push!(metrics, [
            duration,
            num_trips,
            utilization,
            propagated_delay,
            propagated_delay_err,
            distance,
            deadhead,
            geometry
        ])
    end
    
    vehicle_cost = veh_cost * sum(x[1, :])
    link_cost = sum(C .* x) - vehicle_cost
    passenger_cost = sum(propagated_delays) * delay_cost
    passenger_cost_err = sum(propagated_delay_errs) * delay_cost
    service_cost = sum(propagated_delays ./ ridership .* x[2:end, 1]) * op_cost +
        sum(propagated_delays' * (V .* x)[2:end, 2:end]) * op_cost +
        sum(mean(L, dims=2) .* x[:, 1]) * op_cost +
        sum(trips[:, :stop_time] .- trips[:, :start_time]) * op_cost
    service_cost_err = sum(propagated_delay_errs ./ ridership .* x[2:end, 1]) * op_cost +
        sum(propagated_delay_errs' * (V .* x)[2:end, 2:end]) * op_cost

    cost = vehicle_cost + link_cost + passenger_cost + service_cost
    cost_err = passenger_cost_err + service_cost_err
    
    s_10 = ceil(Int, numScenarios/10)
    return SolutionStats(
        sol,
        cost,
        cost_err,
        vehicle_cost,
        link_cost,
        passenger_cost,
        passenger_cost_err,
        service_cost,
        service_cost_err,
        mean(this_s) * 60,
        mean(sort(this_s, dims=2, rev=true)[:, 1:s_10]) * 60,
        std(this_s) * 60,
        std(sort(this_s, dims=2, rev=true)[:, 1:s_10]) * 60,
        1 - total_nis_length / total_duration,
        total_deadhead * 60,
        metrics
    )
end

function getDeadhead(s::Vector{Int}, D::Matrix{Float64})
    time = 0.0

    for i in 1:size(s, 1)-1
        time += D[s[i], s[i+1]]
    end

    return time
end

function getBlockLength(s::Vector{Int}, trips::DataFrame, D::Matrix{Float64})
    start = trips[s[1], :start_time]
    stop = trips[s[end], :stop_time]

    return stop - start + D[1, s[1]] + D[s[end], 1]
end

function notInServiceLength(s::Vector{Int}, trips::DataFrame, D::Matrix{Float64})
    nis = D[1,s[1]]
    nis += D[s[end], 1]

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