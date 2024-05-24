using DataFrames
using Distances
using Distributions
using Graphs
using Random

"""
	VSPInstance

An instance of the Vehicle Scheduling Problem (VSP).

# Fields
- `n::Int`: the number of trips in the network plus one, representing the depot.
- `M::Float64`: maximum possible value for propagated delay of any trip.
- `l::Vector{Distribution}`: normal distribution for expected delays across each trip.
- `C::Matrix{Float64}`: cost for using arc (i, j) in the VSP network.
- `B::Matrix{Float64}`: buffer time on arc (i, j) in the VSP network.
- `D::Matrix{Float64}`: D[i, j] is the distance from the end of trip i to the start of trip j.
- `G::Matrix{Bool}`: 1 if arc (i, j) is usable in the VSP network, 0 otherwise.
- `trips::DataFrame`: the trips DataFrame, see `loadGTFS` in utils.jl.
"""
struct VSPInstance
	n::Int # number of trips + depot
	M::Float64 # big-M constraint variable for propagated delay
	l::Vector{Distribution} # normal distributions for expected delay for all trips
	C::Matrix{Float64} # adjacency-cost lists for all trips
	B::Matrix{Float64} # buffer time between all trips
	D::Matrix{Float64} # D[i, j] is the deadhead time from the end of trip i to the start of trip j
	G::Matrix{Bool} # connections between trips (G[i,j] = 1 if arc i -> j in G)
	trips::DataFrame # trips dataframe
end

"""
	VSPInstance(
		trips::DataFrame[,
		l::Union{Matrix{Float64}, Nothing} = nothing,
		averageSpeed::Float64 = 30.0,
		randomSeed = 1
		]
	)

Create a VSPInstance object from the `trips` DataFrame.

Expected trip delays may be specified by `l` which has mean delays in the first column
and standard deviations in the second.  Possible deadheads are determined
by haversine distance and `averageSpeed`.
"""
function VSPInstance(
	trips::DataFrame;
    l::Union{Matrix{Float64}, Nothing} = nothing, # [μ, σ]
	averageSpeed::Float64 = 30.0,
	randomSeed = 1,
)
	Random.seed!(randomSeed)
	n = size(trips, 1) + 1
	if isnothing(l)
        μ = (trips[:, :stop_time] .- trips[:, :start_time]) * 0.25
        σ = μ * 0.1 / 0.25
        # generate normal distributions for each trip
		l = truncated.(Normal.(μ, σ), 0.0, 4*μ)
	else
		l = truncated.(Normal.(l[:, 1], l[:, 2]), 0.0, 4*μ)
    end
	C = zeros(Float64, n, n)
    C[1, 2:end] .= 100 # cost per vehicle
	B = zeros(Float64, n, n)
	G = zeros(Bool, n, n)
	G[1, 2:end] .= 1 # add link from depot to each trip
	G[2:end, 1] .= 1 # add link from each trip to depot
	D = zeros(Float64, n-1, n-1) # deadhead time between start/stop points

	for i ∈ 1:n-1
		for j ∈ 1:n-1
            i == j && continue
			timeDiff = trips[j, :start_time] - trips[i, :stop_time]
			timeDiff < 0 && continue
			coords_1 = (trips[i, :stop_lat], trips[i, :stop_lon])
			coords_2 = (trips[j, :start_lat], trips[j, :start_lon])
			distance = haversine(coords_1, coords_2, 6372.8)
			D[i, j] = distance / averageSpeed
			if distance / averageSpeed <= timeDiff
				G[i+1, j+1] = 1 # add link if vehicle can deadhead from i -> j
				# ensure rational data
                B[i+1, j+1] = round(timeDiff - distance / averageSpeed; digits = 2)
                C[i+1, j+1] = round(distance + timeDiff; digits = 2)
			end
		end
	end

	# using longest path to tighten big-M
	g = SimpleDiGraph(G[2:end, 2:end])
	m = -1*minimum(Graphs.dijkstra_shortest_paths(g, findall(G[1,2:end]), -1*ones(Int, n-1, n-1)).dists)
	M = sum(sort(maximum.(l), rev=true)[1:m+1])
    B[2:end, 1] .= M # *infinite* buffer time when returning to depot

	return VSPInstance(n, M, l, C, B, D, G, trips)
end

"""
	VSPInstance(inst::VSPInstance, λ::Vector{Float64}, s::Vector{Float64})

Create a VSPInstance object from `inst` with adjusted link costs for Lagrangian
Relaxation.

See `update!` in lagrange.jl.
"""
function VSPInstance(inst::VSPInstance, λ::Vector{Float64}, s::Vector{Float64})
	n = inst.n
	M = inst.M
	l = inst.l
	C = copy(inst.C)
	B = inst.B
	G = inst.G
	trips = inst.trips

	for i ∈ 1:n
		for j ∈ 1:n
			# adjust link costs with respect to s and λ for current iteration
			C[i, j] += λ[j] * (s[i] + l[i] - B[i, j])
		end
	end

	return VSPInstance(n, M, l, C, B, G, trips)
end