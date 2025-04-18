using CSV
using Dates
using DataFrames
using DataFramesMeta
using Distances
using Distributions
using JLD2
using LinearAlgebra
using Random
using Statistics

"""
    string2time(s::String15)

Compute the hour as a Float32 from `s`.

`s` must be formatted as "hh:mm:ss" and string2time supports times beyond "23:59:59".
"""
function string2time(s::String15)
    h, m, s = parse.(Float32, split(s, ":"))
    return h + m / 60 + s / 3600
end

"""
    string2time(dt::Dates.Time)

Compute the hour as a Float32 from the datetime object `dt`.
"""
function string2time(dt::Dates.Time)
    h = hour(dt)
    m = minute(dt)
    s = second(dt)
    return Float32(h + m / 60 + s / 3600)
end

"""
    loadGTFS(path::String)

Load and transform GTFS data located at `path` with respect to the current directory.

Returns two DataFrames.  The first include trip-level information while the second
includes shape-level information.
"""
function loadGTFS(path::String, historical_data)
    # load dataframes for gtfs files
    folder = joinpath(@__DIR__(), path)
    trips = CSV.read(joinpath(folder, "trips.txt"), DataFrames.DataFrame)
    routes = CSV.read(joinpath(folder, "routes.txt"), DataFrames.DataFrame)
    times = CSV.read(joinpath(folder, "stop_times.txt"), DataFrames.DataFrame)
    shapes = CSV.read(joinpath(folder, "shapes.txt"), DataFrames.DataFrame)
    # calendar = CSV.read(joinpath(folder, "calendar_dates.txt"), DataFrames.DataFrame)

    # create trip dataframe ordered by block_id and start_time
    trips_df = @chain trips begin
        # @subset (:service_id .== 2015) # manual filter for Nanaimo
        # @subset (:service_id .== 2086) # manual filter for Cranbrook
        @subset (:service_id .== 595) # manual filter for Victoria
        @rename :direction = :direction_id
        @transform (@byrow :route_id =
            routes[routes.route_id.==:route_id, :].route_short_name[1])
        # @transform(
        #     @byrow :distance = Float32(
        #         last(times[times.trip_id.==:trip_id, :]).shape_dist_traveled / 1000, # convert to km
        #     )
        # )
        @transform (@byrow :start_time =
            string2time(first(times[times.trip_id.==:trip_id, :]).arrival_time))
        @transform (@byrow :stop_time =
            string2time(last(times[times.trip_id.==:trip_id, :]).arrival_time))
        @transform (@byrow :start_lat =
            Float32(first(shapes[shapes.shape_id.==:shape_id, :]).shape_pt_lat))
        @transform (@byrow :start_lon =
            Float32(first(shapes[shapes.shape_id.==:shape_id, :]).shape_pt_lon))
        @transform (@byrow :stop_lat =
            Float32(last(shapes[shapes.shape_id.==:shape_id, :]).shape_pt_lat))
        @transform (@byrow :stop_lon =
            Float32(last(shapes[shapes.shape_id.==:shape_id, :]).shape_pt_lon))
        _[
            :,
            [
                :block_id,
                :route_id,
                :trip_id,
                :shape_id,
                # :distance,
                :direction,
                :start_time,
                :stop_time,
                :start_lat,
                :start_lon,
                :stop_lat,
                :stop_lon,
            ],
        ]
        sort(_, [:block_id, :start_time])
    end

    route_class_dict = Dict(
        row.route => row.route_class for row in eachrow(unique(historical_data, [:route, :route_class]))
    )
    route_dir_dict = Dict(
        row1.route => Dict(
            row2.direction => row2.direction_code for row2 in eachrow(unique(historical_data[historical_data.route .== row1.route, [:direction, :direction_code]]))
        ) for row1 in eachrow(unique(historical_data, [:route]))
    )
    trips_df.route_class = [string(route) in keys(route_class_dict) ? route_class_dict[string(route)] : nothing for route in trips_df.route_id];
    trips_df.direction_code_full = [string(row.route_id) in keys(route_dir_dict) ? 
        row.direction in keys(route_dir_dict[string(row.route_id)]) ? route_dir_dict[string(row.route_id)][row.direction] : nothing : nothing
        for row in eachrow(trips_df)
    ]
    dir_dict = Dict(
        "Inbound" => "IB",
        "Outbound" => "OB",
        "North" => "N",
        "South" => "S",
        "East" => "E",
        "West" => "W",
        "Clockwise" => "CW",
        "Counterclockwise" => "CCW",
        "Counterclo" => "CCW",
        nothing => "?"
    )
    trips_df.direction_code = [dir_dict[code] for code in trips_df.direction_code_full]

    grouped_shapes = groupby(shapes, :shape_id)
    shape_pts = Vector{Vector{Tuple{Float64, Float64}}}()
    shape_dist_traveled = Float64[]

    for shape in grouped_shapes
        # extract shape_pt_lat and shape_pt_lon columns as tuples
        lat_lon_tuples = [(row.shape_pt_lon, row.shape_pt_lat) for row in eachrow(shape)]
        
        # extract shape_dist_traveled
        dist_traveled = Float64[]
        try
            dist_traveled = maximum(shape.shape_dist_traveled) / 1000 # km
        catch
            dist_traveled = sum([haversine(lat_lon_tuples[i], lat_lon_tuples[i+1], 6372.8) for i in 1:length(lat_lon_tuples)-1])
        end
        
        push!(shape_pts, lat_lon_tuples)
        push!(shape_dist_traveled, dist_traveled)
    end

    shapes_df = DataFrame(
        shape_id = unique(shapes.shape_id),
        shape_pts = shape_pts,
        shape_dist_traveled = shape_dist_traveled
    )

    return trips_df, shapes_df
end

"""
    loadHistoricalData(path::String)

Load and historical data located at `path` with respect to the current directory.
"""
function loadHistoricalData(path::String)
    # load dataframe for historical data
    file = joinpath(@__DIR__(), path)
    df = CSV.read(file, DataFrames.DataFrame)

    holidays = [
        Date(2023, 9, 4), # labour day
        Date(2023, 9, 30), # truth and reconciliation
        Date(2023, 10, 9), # thanksgiving
        Date(2023, 11, 11), # remembrance day
        Date(2024, 9, 2), # labour day
        Date(2024, 9, 30), # truth and reconciliation
        Date(2024, 10, 14), # thanksgiving
        Date(2024, 11, 11), # remembrance day
    ]

    is_holiday = [date in holidays for date in df.date]
    if eltype(df.date) != Dates.Date
        df.date = Dates.Date.(df.date, dateformat"m/d/y")
    end
    dow_dict = Dict(
        "Monday" => 1,
        "Tuesday" => 2,
        "Wednesday" => 3,
        "Thursday" => 4,
        "Friday" => 5,
        "Saturday" => 6,
        "Sunday" => 7,
    )
    is_weekday = Dates.dayofweek.(df.date) .> 1 .&& Dates.dayofweek.(df.date) .< 5 .&& [Dates.dayofweek(d) == dow_dict[df.day[i]] for (i, d) in enumerate(df.date)]
    df.dow = Dates.dayofweek.(df.date)
    df.day_from_dict = [dow_dict[d] for d in df.day]
    df = df[.!is_holiday .& is_weekday, :]
    df.hour_of_day = Dates.hour.(df.planned_start_time)
    df.hour_of_day = [h == 0 ? 24 : h for h in df.hour_of_day]
    df.direction_code = df.direction
    df.direction = [(d in ["Inbound", "North", "Clockwise", "East"] ? false : true) for d in df.direction]
    df.planned_start_hour = Dates.hour.(df.planned_start_time) .+ Dates.minute.(df.planned_start_time) ./ 60 .+ Dates.second.(df.planned_start_time) ./ 3600
    df.primary_delay_hours = (df.end_delay_seconds - df.start_delay_seconds) / 3600

    return df
end

"""
    subsetGTFS(
        df::DataFrame[,
        n = nothing,
        routes = nothing,
        start_time = nothing,
        stop_time = nothing,
        randomSeed = 1]
    )

Selects a random subset of trips from `df` of size `n`.

`routes` may be used to select a specific subset of routes.  `start_time` and 
`stop_time` can be used to specify a time frame.
"""
function subsetGTFS(
    df::DataFrame;
    n = nothing,
    routes = nothing,
    start_time = nothing,
    stop_time = nothing,
    randomSeed = 1
)
    subset = copy(df)
    if !isnothing(routes)
        subset = subset[[in(route, routes) for route in subset.route_id], :]
    end
    if !isnothing(start_time)
        subset = subset[[start >= start_time for start in subset.start_time], :]
    end
    if !isnothing(stop_time)
        subset = subset[[stop <= stop_time for stop in subset.stop_time], :]
    end
    if !isnothing(n) && n <= size(subset, 1)
        Random.seed!(randomSeed)
        idx = sample(1:size(subset, 1), n, replace = false)
        subset = subset[idx, :]
    end


    return subset
end

"""
    getDelays(
        trips::DataFrame[,
        randomSeed = nothing,
        meanMulti = 0.1,
        stdMulti = 0.5]
    )


Generates a vector of truncated normal trip delay distributions, one for each trip in `trips`.

The `trips` DataFrame determines how many distributions are generated.  Distribution means are
selected randomly between `(-0.1*trip_lengths, trip_lengths*meanMulti)` and standard deviations
are selected randomly between `(0.1*trip_lengths, trip_lengths.*stdMulti)`.  Distributions are
truncated with an upper bound of `trip_lengths`.
"""
function getDelays(
    trips::DataFrame;
    randomSeed = nothing,
    meanMulti = 0.1,
    stdMulti = 0.5
)
    delays = Distribution[]
    if !isnothing(randomSeed)
        Random.seed!(randomSeed)
    end

    trip_lengths = trips.stop_time - trips.start_time
    μ_dists = Uniform.(-0.1*trip_lengths, trip_lengths*meanMulti)
    σ_dists = Uniform.(0.1*trip_lengths, trip_lengths.*stdMulti)
    μ = rand.(μ_dists)
    σ = rand.(σ_dists)

    for (i, tl) in enumerate(trip_lengths)
  		push!(delays, truncated(Normal(μ[i], σ[i]), upper=tl))
    end

    return delays
end

"""
    getHistoricalDelays(
        trips::DataFrame,
        data::DataFrame,
        n::Int[,
        randomSeed = nothing,
        split = 1.0,
        multi = 1.0]
    )


Generates training and testing delay matrices.

For each trip in `trips`, sample `n` delays from the historical `data`.  `split` determines
the ratio of training to testing samples.  `multi` allows for adjusting delays by a common
factor.  Matches trips to delay data by route number, direction, and trip start hour.
"""
function getHistoricalDelays(
    trips::DataFrame,
    data::DataFrame;
    n = nothing,
    temporal = true,
    randomSeed = nothing,
    split = 0.5,
    new_std = nothing,
    new_mean = nothing,
    new_test_std = nothing,
    new_test_mean = nothing,
    mean_reduc = nothing,
    std_reduc = nothing,
    overlap = 0.0
)
    if !isnothing(randomSeed)
        Random.seed!(randomSeed)
    end
    delays = data.primary_delay_hours
    unique_dates = keys(filter(x -> all(y -> y in x.route, unique(string.(trips.route_id))), groupby(data, :date)))
    unique_dates = [date.date for date in unique_dates]

    n_dates = length(unique_dates)
    if isnothing(n)
        n = n_dates
    end

    if temporal
        unique_dates = sort(unique_dates)
    else
        unique_dates = shuffle(unique_dates)
    end
    n_train_dates = round(Int, n_dates*split)
    n_train = round(Int, n*split)
    n_test = n-n_train
    n_offset = round(Int, overlap*n_test)

    if temporal
        train_dates = unique_dates[year.(unique_dates) .== 2023]
    elseif n_train > n_train_dates
        train_dates = sample(unique_dates[1:n_train_dates], n_train, replace=true)
    else
        train_dates = sample(unique_dates[1:n_train_dates], n_train, replace=false)
    end
    L_train = zeros(Float64, size(trips, 1), length(train_dates))

    if temporal
        test_dates = unique_dates[year.(unique_dates) .== 2024]
        n_test = length(test_dates)
        n_offset = round(Int, overlap*n_test)
        test_dates = sample(test_dates, n_test-n_offset, replace=false)
        test_dates = append!(test_dates, sample(train_dates, n_offset, replace=false))
    elseif n_test > n_dates-n_train_dates
        test_dates = sample(unique_dates[n_train_dates+1-n_offset:end-n_offset], n_test, replace=true)
    else
        test_dates = sample(unique_dates[n_train_dates+1-n_offset:end-n_offset], n_test, replace=false)
    end
    L_test = zeros(Float64, size(trips, 1), length(test_dates))

    for (i, trip) in enumerate(eachrow(trips))
        for (j, date) in enumerate(train_dates)
            route_str = string(trip.route_id)
            mask = (data.route .== route_str) .& (data.direction .== trip.direction) .& (data.date .== date)

            if sum(mask) < 1
                mask = (data.route .== route_str) .& (data.date .== date)
            end
            if sum(mask) < 1
                error("Route $(route_str) at $(trip.start_time) returned no training results for $date.")
            end

            subset = delays[mask]
            start_times = data.planned_start_hour[mask]
            closest_index = argmin(abs.(start_times .- trip.start_time))

            L_train[i, j] = subset[closest_index]
        end

        for (j, date) in enumerate(test_dates)
            route_str = string(trip.route_id)
            mask = (data.route .== route_str) .& (data.direction .== trip.direction) .& (data.date .== date)

            if sum(mask) < 1
                mask = (data.route .== route_str) .& (data.date .== date)
            end
            if sum(mask) < 1
                error("Route $(route_str) at $(trip.start_time) returned no training results for $date.")
            end

            subset = delays[mask]
            start_times = data.planned_start_hour[mask]
            closest_index = argmin(abs.(start_times .- trip.start_time))

            L_test[i, j] = subset[closest_index]
        end
    end

    train_mean = mean(L_train)
    train_std = std(L_train)
    test_mean = mean(L_test)
    test_std = std(L_test)
    
    if !isnothing(new_std)
        L_train = (L_train .- train_mean) ./ train_std .* new_std ./ 60 .+ train_mean
        L_test = (L_test .- test_mean) ./ test_std .* new_std ./ 60 .+ test_mean
    elseif !isnothing(new_test_std)
        L_test = (L_test .- test_mean) ./ test_std .* new_test_std ./ 60 .+ test_mean
    elseif !isnothing(std_reduc)
        L_train = (L_train .- train_mean) .* (1-std_reduc) .+ train_mean
        L_test = (L_test .- test_mean) .* (1-std_reduc) .+ test_mean
    end

    if !isnothing(new_mean)
        L_train = (L_train .- train_mean) .+ (new_mean / 60)
        L_test = (L_test .- test_mean) .+ (new_mean / 60)
    elseif !isnothing(new_test_mean)
        L_test = (L_test .- test_mean) .+ (new_test_mean / 60)
    elseif !isnothing(mean_reduc)
        L_train[L_train .> 0] .= L_train[L_train .> 0] .* (1-mean_reduc)
        L_test[L_test .> 0] .= L_test[L_test .> 0] .* (1-mean_reduc)
    end

    return L_train, L_test
end

"""
    getRidership(
        trips::DataFrame[,
        randomSeed = nothing,
        min = 0,
        max = 120,]
    )

Generates a vector representing average ridership across each trip in `trips`.

Each ridership value is selected randomly between `min` and `max`.
"""
function getRidership(
    trips::DataFrame;
    randomSeed = nothing,
    min = 0,
    max = 120,
)
    if !isnothing(randomSeed)
        Random.seed!(randomSeed)
    end

    n = size(trips, 1)
    dist = Uniform(min, max)

    return rand(dist, n)
end

"""
    getHistoricalRidership(
        trips::DataFrame,
        data::DataFrame
    )

Generates a vector representing average ridership across each trip in `trips`.

Each ridership value is the mean of historical `data`.
Matches trips to delay data by route number, direction, and trip start hour.
"""
function getHistoricalRidership(
    trips::DataFrame,
    data::DataFrame
)
    r = Float64[]

    for trip in eachrow(trips)
        route_str = string(trip.route_id)
        subset = data.total_boardings[
            (data.route .== route_str) .&
            (abs.(data.planned_start_hour .- trip.start_time) .< 1) .&
            (data.direction .== trip.direction)
        ]

        if isempty(subset)
            subset = data.total_boardings[
                (data.route .== route_str) .&
                (abs.(data.planned_start_hour .- trip.start_time) .< 2) .&
                (data.direction .== trip.direction)
            ]
        end

        push!(r, max(mean(subset), 1))
    end

    return r
end

# Function to check if two line segments (p1, q1) and (p2, q2) intersect
function do_segments_intersect(p1, q1, p2, q2)
    # Helper function to find the orientation of the ordered triplet (p, q, r)
    # The function returns:
    # 0 -> p, q and r are collinear
    # 1 -> Clockwise
    # 2 -> Counterclockwise
    function orientation(p, q, r)
        val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
        if val == 0
            return 0  # collinear
        elseif val > 0
            return 1  # clockwise
        else
            return 2  # counterclockwise
        end
    end

    # Helper function to check if point q lies on line segment pr
    function on_segment(p, q, r)
        return q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]) &&
               q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2])
    end

    # Find the four orientations needed for the general and special cases
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if o1 != o2 && o3 != o4
        return true
    end

    # Special cases
    # p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if o1 == 0 && on_segment(p1, p2, q1)
        return true
    end

    # p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if o2 == 0 && on_segment(p1, q2, q1)
        return true
    end

    # p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if o3 == 0 && on_segment(p2, p1, q2)
        return true
    end

    # p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if o4 == 0 && on_segment(p2, q1, q2)
        return true
    end

    # Doesn't fall in any of the above cases
    return false
end

# Main function to check if two vectors of lat/lon tuples intersect
function do_paths_intersect(path1::Vector{Tuple{Float64, Float64}}, path2::Vector{Tuple{Float64, Float64}})
    n1 = length(path1)
    n2 = length(path2)
    
    for i in 1:(n1-1)
        for j in 1:(n2-1)
            if do_segments_intersect(path1[i], path1[i+1], path2[j], path2[j+1])
                return true
            end
        end
    end
    return false
end