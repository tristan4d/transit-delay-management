using CSV
using Dates
using DataFrames
using DataFramesMeta
using Distributions

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
function loadGTFS(path::String)
    # load dataframes for gtfs files
    folder = joinpath(@__DIR__(), path)
    trips = CSV.read(joinpath(folder, "trips.txt"), DataFrames.DataFrame)
    routes = CSV.read(joinpath(folder, "routes.txt"), DataFrames.DataFrame)
    times = CSV.read(joinpath(folder, "stop_times.txt"), DataFrames.DataFrame)
    shapes = CSV.read(joinpath(folder, "shapes.txt"), DataFrames.DataFrame)
    # calendar = CSV.read(joinpath(folder, "calendar_dates.txt"), DataFrames.DataFrame)

    # create trip dataframe ordered by block_id and start_time
    trips_df = @chain trips begin
        @subset (:service_id .== 2015) # manual filter for Nanaimo
        # @subset (:service_id .== 2086) # manual filter for Cranbrook
        @transform (@byrow :route_id =
            routes[routes.route_id.==:route_id, :].route_short_name[1])
        @transform(
            @byrow :distance = Float32(
                last(times[times.trip_id.==:trip_id, :]).shape_dist_traveled / 1000, # convert to km
            )
        )
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
                :distance,
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

    grouped_shapes = groupby(shapes, :shape_id)
    shape_pts = Vector{Vector{Tuple{Float64, Float64}}}()
    shape_dist_traveled = Float64[]

    for shape in grouped_shapes
        # extract shape_pt_lat and shape_pt_lon columns as tuples
        lat_lon_tuples = [(row.shape_pt_lon, row.shape_pt_lat) for row in eachrow(shape)]
        
        # extract shape_dist_traveled
        dist_traveled = sum(shape.shape_dist_traveled)
        
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
    subset = df
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