using CSV
using Dates
using DataFrames
using DataFramesMeta
using Distributions
using LinearAlgebra

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
        dist_traveled = maximum(shape.shape_dist_traveled) / 1000 # km
        
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

function primaryDelays(trips::DataFrame; form = 0, bbox = nothing, shapes = nothing)
    n = size(trips, 1)
    min_start = minimum(trips.start_time)
    max_start = maximum(trips.start_time)
    l = zeros(Float64, n, 2)
    t = trips[:, :stop_time] - trips[:, :start_time]

    if form == 0
        dist = Uniform(min_start, max_start)
        x = trips.start_time
    elseif form == 1
        dist = Chi(1)
        x = (trips.start_time .- min_start) ./ (max_start - min_start) .* 2
    elseif form == 2
        dist = Chi(1)
        x = ((trips.start_time .- min_start) ./ (max_start - min_start) .- 1) .* (-2)
    end

    λ = pdf.(dist, x) / 2
    l[:, 1] = t.*λ
    l[:, 2] = t.*λ/2

    if !isnothing(bbox) && !isnothing(shapes)
        for i in 1:n
            if do_paths_intersect(shapes[shapes.shape_id .== trips[i, :shape_id], :shape_pts][1], bbox)
                l[i, 1] *= 2
            end
        end
    end

    return l
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