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
    extra_cost::Float64
    secondary_delays::Matrix{Float64}
    end_of_trip_delays::Matrix{Float64}
    utilization::Float64
    deadhead::Float64
    buffer::Float64
    buffer_std::Float64
    buffer_err::Tuple{Float64, Float64}
    tt::Float64
    tt_std::Float64
    tt_err::Tuple{Float64, Float64}
    # metrics::DataFrame
end

function getSolutionStats(
    sol::VSPSolution;
    Δ = nothing,
    delays = nothing,
    ridership = nothing,
    extra_costs = 0.0,
    base_vehicles = nothing
)
    this_x = convert(Matrix{Bool}, round.(sol.x))
    if isnothing(Δ)
        Δ = sol.Δ
    end
    schedules = generate_blocks(this_x)
    n = sol.mod.inst.n
    op_cost = sol.mod.inst.op_cost
    idle_cost = sol.mod.inst.idle_cost
    delay_cost = sol.mod.inst.delay_cost
    veh_cost = sol.mod.inst.veh_cost
    tts_ratio = sol.mod.inst.tts_ratio
    rta_mask = sol.mod.inst.rta_mask
    B = sol.mod.inst.B
    Q = sol.mod.inst.Q
    C = sol.mod.inst.C .- Q
    D = sol.mod.inst.D
    V = sol.mod.inst.V
    trips = sol.mod.inst.trips

    if isnothing(delays)
        L = sol.mod.inst.L_test
    else
        L = delays
    end
    if isnothing(ridership)
        ridership = sol.mod.inst.r
    end

    this_s = sol.s .* ridership

    total_duration = 0.0
    total_nis_length = 0.0
    total_deadhead = 0.0
    for schedule in schedules
        duration = getBlockLength(schedule, Δ, trips, D, V)
        total_duration += duration
        nis_length = notInServiceLength(schedule, Δ, trips, D, V)
        total_nis_length += nis_length
        deadhead = getDeadhead(schedule, D)
        total_deadhead += deadhead
    end
    buffer_mask = .!(V[2:end, 2:end])
    true_buffers = vec((B[2:end, 2:end] .- Δ)[this_x[2:end, 2:end] .* buffer_mask]) .* 60
    
    vehicle_cost = veh_cost * sum(this_x[1, :])
    link_cost = sum(C .* this_x) - vehicle_cost + sum(Δ .* op_cost)
    vehicle_cost -= veh_cost * (isnothing(base_vehicles) ? 0 : base_vehicles)
    # service_cost = sum(sum(this_x[2:end, :], dims=2) .* (trips[:, :stop_time] .- trips[:, :start_time])) * op_cost
    service_cost = sum(
            trips[1:length(rta_mask), :stop_time] .- trips[1:length(rta_mask), :start_time] .+ max.(
                Δ, mean(max.(L, 0); dims=1)
            )
        ) * op_cost
    passenger_tt_cost = sum(m .* tts_ratio .* delay_cost .* ridership)
    # passenger_cost = sum(max.(propagated_delays, 0)) * delay_cost
    passenger_cost = sum(this_s) * delay_cost
    passenger_cost_std = sum(sol.s_std .* ridership) * delay_cost
    # passenger_cost_lo = sum(max.(propagated_delay_lo, 0)) * delay_cost
    passenger_cost_lo = sum(sol.s_lo .* ridership) * delay_cost
    # passenger_cost_hi = sum(max.(propagated_delay_hi, 0)) * delay_cost
    passenger_cost_hi = sum(sol.s_hi .* ridership) * delay_cost

    cost = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost + extra_costs
    cost_lo = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost_lo + extra_costs
    cost_hi = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost_hi + extra_costs
    
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
        extra_costs,
        B,
        B,
        1 - total_nis_length / total_duration,
        total_deadhead * 60,
        mean(true_buffers),
        std(true_buffers),
        (quantile(true_buffers, 0.25), quantile(true_buffers, 0.75)),
        mean(Δ),
        std(Δ),
        (quantile(Δ, 0.25), quantile(Δ, 0.75))
    )
end

function getSolutionStats(
    x::Union{Matrix{Float64}, Matrix{Bool}},
    instance::VSPInstance;
    m = nothing,
    delays = nothing,
    ridership = nothing,
    extra_costs = 0.0,
    base_vehicles = nothing
)
    this_x = convert(Matrix{Bool}, round.(x))
    if isnothing(m)
        m = zeros(Float64, instance.n-1)
    end
    schedules = generate_blocks(this_x)
    n = instance.n
    op_cost = instance.op_cost
    idle_cost = instance.idle_cost
    delay_cost = instance.delay_cost
    veh_cost = instance.veh_cost
    tts_ratio = instance.tts_ratio
    rta_mask = instance.rta_mask
    B = instance.B
    Q = instance.Q
    C = instance.C .- Q
    D = instance.D
    V = instance.V
    trips = instance.trips

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
    secondary_delays = zeros(Float64, n, numScenarios)
    end_of_trip_delays = zeros(Float64, n, numScenarios)

    for scenario in 1:numScenarios
        secondary_delays[:, scenario] = feasibleDelays(this_x, L[:, scenario], B, m=m)
        end_of_trip_delays[:, scenario] = feasibleDelays(this_x, L[:, scenario], B, m=m)
        this_s[2:end, scenario] .= end_of_trip_delays[2:end, scenario] .* ridership
    end

    total_duration = 0.0
    total_nis_length = 0.0
    total_deadhead = 0.0
    for schedule in schedules
        duration = getBlockLength(schedule, m, trips, D, V)
        total_duration += duration
        nis_length = notInServiceLength(schedule, m, trips, D, V)
        total_nis_length += nis_length
        deadhead = getDeadhead(schedule, D)
        total_deadhead += deadhead
    end
    buffer_mask = .!(V[2:end, 2:end])
    true_buffers = vec((B[2:end, 2:end] .- m)[this_x[2:end, 2:end] .* buffer_mask]) .* 60
    
    vehicle_cost = veh_cost * sum(this_x[1, :])
    link_cost = sum(C .* this_x) - vehicle_cost + sum(m .* op_cost)
    vehicle_cost -= veh_cost * (isnothing(base_vehicles) ? 0 : base_vehicles)
    # service_cost = sum(sum(this_x[2:end, :], dims=2) .* (trips[:, :stop_time] .- trips[:, :start_time])) * op_cost
    service_cost = sum(
            trips[1:length(rta_mask), :stop_time] .- trips[1:length(rta_mask), :start_time] .+ max.(
                m, mean(max.(L, 0); dims=1)
            )
        ) * op_cost
    passenger_tt_cost = sum(m .* tts_ratio .* delay_cost .* ridership)
    # passenger_cost = sum(max.(propagated_delays, 0)) * delay_cost
    passenger_cost = mean(sum(this_s, dims=1)) * delay_cost
    passenger_cost_std = std(sum(this_s, dims=1)) * delay_cost
    # passenger_cost_lo = sum(max.(propagated_delay_lo, 0)) * delay_cost
    passenger_cost_lo = quantile(vec(sum(this_s, dims=1)), 0.25) * delay_cost
    # passenger_cost_hi = sum(max.(propagated_delay_hi, 0)) * delay_cost
    passenger_cost_hi = quantile(vec(sum(this_s, dims=1)), 0.75) * delay_cost

    cost = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost + extra_costs
    cost_lo = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost_lo + extra_costs
    cost_hi = vehicle_cost + link_cost + passenger_tt_cost + passenger_cost_hi + extra_costs
    
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
        extra_costs,
        secondary_delays,
        end_of_trip_delays,
        1 - total_nis_length / total_duration,
        total_deadhead * 60,
        mean(true_buffers),
        std(true_buffers),
        (quantile(true_buffers, 0.25), quantile(true_buffers, 0.75)),
        mean(m),
        std(m),
        (quantile(m, 0.25), quantile(m, 0.75))
    )
end

function getDeadhead(s::Vector{Int}, D::Matrix{Float64})
    time = 0.0

    for i in 1:size(s, 1)-1
        time += D[s[i]+1, s[i+1]+1]
    end

    return time
end

function getBlockLength(s::Vector{Int}, m::Vector{Float64}, trips::DataFrame, D::Matrix{Float64}, V::Matrix{Bool})
    start = trips[s[1], :start_time]
    stop = trips[s[end], :stop_time] + m[s[end]]

    total_length = stop - start + D[1, s[1]+1] + D[s[end]+1, 1]

    for i in 1:length(s)-1
        if V[s[i]+1, s[i+1]+1]
            start = trips[s[i], :stop_time] + m[s[i]]
            stop = trips[s[i+1], :start_time]
            total_length -= stop - start + D[s[i]+1, s[i+1]+1]
        end
    end

    return total_length
end

function notInServiceLength(s::Vector{Int}, m::Vector{Float64}, trips::DataFrame, D::Matrix{Float64}, V::Matrix{Bool})
    nis = D[1,s[1]+1]
    nis += D[s[end]+1, 1]

    for i in 1:length(s)-1
        V[s[i]+1, s[i+1]+1] && continue
        start = trips[s[i], :stop_time] + m[s[i]]
        stop = trips[s[i+1], :start_time]

        nis += stop - start
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