using DataFrames

struct SolutionStats
    sol::Union{VSPSolution, MCFSolution}
    metrics::DataFrame
end

function SolutionStats(sol::Union{VSPSolution, MCFSolution})
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
            Float64[]
        ],
        [
            "duration",
            "num_trips",
            "utilization",
            "propagated_delay"
        ]
    )

    for schedule in schedules
        duration = getBlockLength(schedule, trips)
        num_trips = length(schedule)
        utilization = 1 - notInServiceLength(schedule, trips) / duration
        propagated_delay = sum(propagated_delays[schedule])
        push!(metrics, [duration, num_trips, utilization, propagated_delay])
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