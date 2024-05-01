using CSV
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
    loadGTFS(path::String)

Load and transform GTFS data located at `path` with respect to the current directory.
"""
function loadGTFS(path::String)
    # load dataframes for gtfs files
    folder = joinpath(@__DIR__(), path)
    trips = CSV.read(joinpath(folder, "trips.txt"), DataFrames.DataFrame)
    routes = CSV.read(joinpath(folder, "routes.txt"), DataFrames.DataFrame)
    times = CSV.read(joinpath(folder, "stop_times.txt"), DataFrames.DataFrame)
    shapes = CSV.read(joinpath(folder, "shapes.txt"), DataFrames.DataFrame)
    calendar = CSV.read(joinpath(folder, "calendar_dates.txt"), DataFrames.DataFrame)

    # create trip dataframe ordered by block_id and start_time
    df = @chain trips begin
        @subset (:service_id .== 2015) # manual filter for Tuesday
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

    return df
end

"""
    subsetGTFS(df::DataFrame, n::Int[, randomSeed = 1])

Selects a random subset of trips from `df` of size `n`.
"""
function subsetGTFS(df::DataFrame, n::Int; randomSeed = 1)
    Random.seed!(randomSeed)
    idx = sample(1:size(df, 1), n, replace = false)

    return df[idx, :]
end