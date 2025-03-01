using DataFrames
using GeoInterface
using Statistics

"""
    SolutionStats

Solution statistics and metrics for minimum cost flow or delay-aware models.

# Fields
- `cost::Float64`: cost, in monetary units, of the solution.
- `cost_err::Tuple{Float64}`: cost confidence interval, in monetary units, of the solution.
- `vehicle_cost::Float64`: cost, in monetary units, of the fleet.
- `link_cost::Float64`: cost, in monetary units, of the chosen links.
- `service_cost::Float64`: cost, in monetary units, of the service delivered.
- `passenger_cost::Float64`: cost, in monetary units, of the total passenger delay.
- `passenger_cost_err::Tuple{Float64}`: cost confidence interval, in monetary units, of the total passenger delay.
- `utilization::Float64`: the percent of time spent moving passengers.
- `deadhead::Float64`: the distance, in minutes, of deadheading in the solution.
- `metrics::DataFrame`: trip-level metrics.
"""
struct SolutionStats
    cost::Float64
    cost_err::Tuple{Float64, Float64}
    vehicle_cost::Float64
    link_cost::Float64
    service_cost::Float64
    passenger_tt_cost::Float64
    passenger_cost::Float64
    passenger_cost_std::Float64
    passenger_cost_err::Tuple{Float64, Float64}
    utilization::Float64
    deadhead::Float64
    metrics::DataFrame
end

function getSolutionStats(
    x::Union{Matrix{Float64}, Matrix{Bool}},
    instance::VSPInstance;
    shapes = nothing,
    endoftrip = true,
    delays = nothing,
    ridership = nothing
)
    this_x = convert(Matrix{Bool}, round.(x))
    schedules = generate_blocks(this_x)
    n = instance.n
    op_cost = instance.op_cost
    delay_cost = instance.delay_cost
    veh_cost = instance.veh_cost
    rta_mask = instance.rta_mask
    B = instance.B
    Q = instance.Q
    C = instance.C .- Q
    D = instance.D
    trips = instance.trips
    propagated_delays = zeros(Float64, n-1)

    if isnothing(delays)
        L = instance.L_test
    else
        L = delays
    end
    if isnothing(ridership)
        ridership = instance.r
    end

    numScenarios = size(L, 2)
    this_s = zeros(Float64, n, numScenarios)

    for scenario in 1:numScenarios
        this_s[:, scenario] = feasibleDelays(this_x, L[:, scenario], B, endoftrip=endoftrip)
        this_s[2:end, scenario] .*= ridership
    end

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

    propagated_delays = vec(mean(this_s, dims=2))[2:end]
    propagated_delay_lo = vec(quantile.(eachrow(this_s), .25))[2:end]
    propagated_delay_hi = vec(quantile.(eachrow(this_s), .75))[2:end]
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
    
    vehicle_cost = veh_cost * sum(this_x[1, :])
    link_cost = sum(C .* this_x) - vehicle_cost
    # service_cost = sum(sum(this_x[2:end, :], dims=2) .* (trips[:, :stop_time] .- trips[:, :start_time])) * op_cost
    service_cost = sum(trips[1:length(rta_mask), :stop_time] .- trips[1:length(rta_mask), :start_time]) * op_cost
    passenger_tt_cost = sum(Q .* this_x)
    # passenger_cost = sum(max.(propagated_delays, 0)) * delay_cost
    passenger_cost = mean(sum(this_s, dims=1)) * delay_cost
    passenger_cost_std = std(sum(this_s, dims=1)) * delay_cost
    # passenger_cost_lo = sum(max.(propagated_delay_lo, 0)) * delay_cost
    passenger_cost_lo = quantile(vec(sum(this_s, dims=1)), 0.25) * delay_cost
    # passenger_cost_hi = sum(max.(propagated_delay_hi, 0)) * delay_cost
    passenger_cost_hi = quantile(vec(sum(this_s, dims=1)), 0.75) * delay_cost

    cost = vehicle_cost + link_cost + service_cost + passenger_tt_cost + passenger_cost
    cost_lo = vehicle_cost + link_cost + service_cost + passenger_tt_cost + passenger_cost_lo
    cost_hi = vehicle_cost + link_cost + service_cost + passenger_tt_cost + passenger_cost_hi
    
    return SolutionStats(
        cost,
        (cost_lo, cost_hi),
        vehicle_cost,
        link_cost,
        service_cost,
        passenger_tt_cost,
        passenger_cost,
        passenger_cost_std,
        (passenger_cost_lo, passenger_cost_hi),
        1 - total_nis_length / total_duration,
        total_deadhead,
        metrics
    )
end

function getDeadhead(s::Vector{Int}, D::Matrix{Float64})
    time = 0.0

    for i in 1:size(s, 1)-1
        time += D[s[i]+1, s[i+1]+1]
    end

    return time
end

function getBlockLength(s::Vector{Int}, trips::DataFrame, D::Matrix{Float64})
    start = trips[s[1], :start_time]
    stop = trips[s[end], :stop_time]

    return stop - start + D[1, s[1]+1] + D[s[end]+1, 1]
end

function notInServiceLength(s::Vector{Int}, trips::DataFrame, D::Matrix{Float64})
    nis = D[1,s[1]+1]
    nis += D[s[end]+1, 1]

    for i in 1:length(s)-1
        stop = trips[s[i], :stop_time]
        start = trips[s[i+1], :start_time]

        nis += start - stop
    end

    return nis
end

function getGeometry(
    s::Vector{Int},
    trips::DataFrame,
    shapes::Union{DataFrame, Nothing}
)
    coordinates = GeoInterface.LineString[]
    distance = 0.0
    if isnothing(shapes)
        return distance, coordinates
    end
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
    xs::Vector{Matrix{Bool}}
)
    return sum(reduce(.&, xs)) / sum(reduce(.|, xs))
end