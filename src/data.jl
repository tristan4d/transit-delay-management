using DataFrames
using Distances
using Distributions
using Graphs

"""
	VSPInstance

An instance of the Vehicle Scheduling Problem (VSP).

# Fields
- `n::Int`: the number of trips in the network plus one, representing the depot.
- `M::Float64`: maximum possible value for propagated delay of any trip.
- `op_cost::Int`: the cost per hour of service operation.
- `delay_cost::Int`: the cost per hour that a passenger is delayed.
- `veh_cost::Int`: the daily cost of adding a new vehicle to the solution.
- `l::Vector{Distribution}`: normal distribution for expected delays across each trip.
- `r::Vector{Float64}`: ridership for each trip.
- `C::Matrix{Float64}`: cost for using arc (i, j) in the VSP network.
- `B::Matrix{Float64}`: buffer time on arc (i, j) in the VSP network.
- `D::Matrix{Float64}`: D[i, j] is the distance from the end of trip i to the start of trip j.
- `V::Matrix{Bool}`: 1 if arc (i, j) returns to the depot in between.
- `G::Matrix{Bool}`: 1 if arc (i, j) is usable in the VSP network, 0 otherwise.
- `trips::DataFrame`: the trips DataFrame, see `loadGTFS` in utils.jl.
"""
struct VSPInstance
	n::Int # number of trips + depot
	M::Float64 # big-M constraint variable for propagated delay
	op_cost::Int # $ per hour of operation
	delay_cost::Int # $ per passenger waiting hour
	veh_cost::Int # $ per vehicle
	l::Vector{Distribution} # normal distributions for expected delay for all trips
	r::Vector{Float64} # ridership for all trips
	C::Matrix{Float64} # adjacency-cost lists for all trips
	B::Matrix{Float64} # buffer time between all trips
	D::Matrix{Float64} # D[i, j] is the deadhead time from the end of trip i to the start of trip j
	V::Matrix{Bool} # virtual depot trips
	G::Matrix{Bool} # connections between trips (G[i,j] = 1 if arc i -> j in G)
	trips::DataFrame # trips dataframe
end

"""
	VSPInstance(
		trips::DataFrame[,
		randomSeed = 1,
		l::Union{Matrix{Float64}, Vector{Distribution}, Nothing} = nothing,
		r::Union{Vector{Float64}, Nothing} = nothing,
		op_cost = 160,
		delay_cost = 37,
		veh_cost = 600,
		averageSpeed::Float64 = 30.0,
		depot_loc = (mean(trips[:, :start_lat]), mean(trips[:, :start_lon]))
		]
	)

Create a VSPInstance object from the `trips` DataFrame.

Expected trip delays may be specified by `l` which is a vector of distributions or 
a matrix with mean delays in the first columnand standard deviations in the second.
Trip ridership may be specified by `r`.  The costing of solutions may be modified with

- `op_cost`, the cost of one hour of service operations;
- `delay_cost`, the cost of an hour of passenger delay; and
- `veh_cost`, the daily cost of adding a new vehicle to the solution.

Possible deadheads are determined by haversine distance and `averageSpeed`.  A `depot_loc`
may be specified, and defaults to the average starting location of all trips in `trips`.
"""
function VSPInstance(
	trips::DataFrame;
	randomSeed = 1,
    l::Union{Matrix{Float64}, Vector{Distribution}, Nothing} = nothing, # [μ, σ]
	r::Union{Vector{Float64}, Nothing} = nothing,
	op_cost = 160, # $ per hour of operation
	delay_cost = 37, # $ per passenger waiting hour
	veh_cost = 600, # $ per vehicle
	averageSpeed::Float64 = 30.0,
	depot_loc = (mean(trips[:, :start_lat]), mean(trips[:, :start_lon]))
)
	n = size(trips, 1) + 1
	delays = Distribution[]
	if isnothing(l)
		delays = getDelays(trips; randomSeed=randomSeed)
	else
		if eltype(l) == Distribution
			delays = l
		else
			trip_lengths = trips[:, :stop_time] .- trips[:, :start_time]
			μ = l[:, 1]
			σ = l[:, 2]

			for (i, tl) in enumerate(trip_lengths)
				push!(delays, truncated(Normal(μ[i], σ[i]), upper=tl))
			end
		end
	end
	C = zeros(Float64, n, n)
    C[1, 2:end] .= veh_cost # cost per vehicle
	B = zeros(Float64, n, n)
	G = zeros(Bool, n, n)
	G[1, 2:end] .= 1 # add link from depot to each trip
	G[2:end, 1] .= 1 # add link from each trip to depot
	D = zeros(Float64, n-1, n-1) # deadhead time between start/stop points
	V = zeros(Bool, n, n)

	for i ∈ 1:n-1
		D[1, i] = haversine(depot_loc, (trips[i, :start_lat], trips[i, :start_lon]), 6372.8) / averageSpeed 
		D[i, 1] = haversine((trips[i, :stop_lat], trips[i, :stop_lon]), depot_loc, 6372.8) / averageSpeed 
		C[1, i] += D[1, i] * op_cost
		C[i, 1] = D[i, 1] * op_cost
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
				if timeDiff < 3.0
	                C[i+1, j+1] = round(distance / averageSpeed + timeDiff; digits = 2) * op_cost # cost per hour
				else
					V[i+1, j+1] = 1
					d1 = haversine(coords_1, depot_loc, 6372.8)
					d2 = haversine(depot_loc, coords_2, 6372.8)
					D[i, j] = (d1+d2) / averageSpeed
					C[i+1, j+1] = round((d1 + d2) / averageSpeed; digits = 2) * op_cost * 2 # return to depot
				end
			end
		end
	end

	# using longest path to tighten big-M
	g = SimpleDiGraph(G[2:end, 2:end])
	m = -1*minimum(Graphs.dijkstra_shortest_paths(g, findall(G[1,2:end]), -1*ones(Int, n-1, n-1)).dists)
	M = sum(sort(maximum.(delays), rev=true)[1:m+1])
    B[2:end, 1] .= M # *infinite* buffer time when returning to depot

	# get ridership
	if isnothing(r)
		r = getRidership(trips; randomSeed=randomSeed)
	end

	return VSPInstance(n, M, op_cost, delay_cost, veh_cost, delays, r, C, B, D, V, G, trips)
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