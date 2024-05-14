using DataFrames
using GeoInterface

struct SolutionStats
    sol::Union{VSPSolution, MCFSolution}
    metrics::DataFrame
end

function getSolutionStats(sol::Union{VSPSolution, MCFSolution}, shapes::DataFrame)
    x = convert(Matrix{Bool}, round.(sol.x))
    schedules = generate_blocks(x)
    l = sol.mod.inst.l
    B = sol.mod.inst.B
    propagated_delays = feasibleDelays(sol.x, l, B)[2:end]
    trips = sol.mod.inst.trips
    metrics = DataFrame(
        [
            Float64[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Vector{Vector{GeoInterface.LineString}}()
        ],
        [
            "duration",
            "num_trips",
            "utilization",
            "propagated_delay",
            "trip_distance",
            "geometry"
        ]
    )

    for schedule in schedules
        duration = getBlockLength(schedule, trips)
        num_trips = length(schedule)
        utilization = 1 - notInServiceLength(schedule, trips) / duration
        propagated_delay = sum(propagated_delays[schedule])
        distance, geometry = getGeometry(schedule, trips, shapes)
        push!(metrics, [
            duration,
            num_trips,
            utilization,
            propagated_delay,
            distance,
            geometry
        ])
    end
    
    return SolutionStats(sol, metrics)
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
        try
            dh_start = GeoInterface.getpoint(coordinates[t-1], GeoInterface.npoint(coordinates[t-1]))
            dh_end = GeoInterface.getpoint(geom, 1)
            push!(coordinates, GeoInterface.LineString([dh_start, dh_end]))
        catch
            nothing
        end

        push!(coordinates, geom)
        distance += shapes[shapes.shape_id .== trips[trip, :shape_id], :shape_dist_traveled][1]
    end

    return distance, coordinates
end