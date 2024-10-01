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
- `L_train::Matrix{Float64}`: delay matrix for training.
- `L_test::Matrix{Float64}`: delay matrix for testing.
- `r::Vector{Float64}`: ridership for each trip.
- `C::Matrix{Float64}`: cost for using arc (i, j) in the VSP network.
- `B::Matrix{Float64}`: buffer time on arc (i, j) in the VSP network.
- `D::Matrix{Float64}`: D[i, j] is the distance from the end of trip i to the start of trip j.
- `V::Matrix{Bool}`: 1 if arc (i, j) returns to the depot in between.
- `G::Matrix{Bool}`: 1 if arc (i, j) is usable in the VSP network, 0 otherwise.
- `trips::DataFrame`: the trips DataFrame, see `loadGTFS` in utils.jl.
"""
struct VSPInstance
	n::Int
	M::Float64
	op_cost::Int
	delay_cost::Int
	veh_cost::Int
	L_train::Matrix{Float64}
	L_test::Matrix{Float64}
	r::Vector{Float64}
	C::Matrix{Float64}
	B::Matrix{Float64}
	D::Matrix{Float64}
	V::Matrix{Bool}
	G::Matrix{Bool}
	trips::DataFrame
end

"""
	VSPInstance(
		trips::DataFrame,
		randomSeed = 1,
		L_train::Matrix{Float64},
		L_test::Matrix{Float64}[,
		op_cost = 160,
		delay_cost = 37,
		veh_cost = 600,
		averageSpeed::Float64 = 30.0,
		depot_loc = (mean(trips[:, :start_lat]), mean(trips[:, :start_lon]))]
	)

Create a VSPInstance object from the `trips` DataFrame.

Trip ridership may be specified by `r`.
`L_train` and `L_test` are the training and testing delays for this instance.
The costing of solutions may be modified with

- `op_cost`, the cost of one hour of service operations;
- `delay_cost`, the cost of an hour of passenger delay; and
- `veh_cost`, the daily cost of adding a new vehicle to the solution.

Possible deadheads are determined by haversine distance and `averageSpeed`.  A `depot_loc`
may be specified, and defaults to the average starting location of all trips in `trips`.
"""
function VSPInstance(
	trips::DataFrame,
	r::Vector{Float64},
	L_train::Matrix{Float64},
	L_test::Matrix{Float64};
	original = false,
	rta_only = false,
	op_cost = 160, # $ per hour of operation
	delay_cost = 37, # $ per passenger waiting hour
	veh_cost = 600, # $ per vehicle
	averageSpeed::Float64 = 30.0,
	depot_loc = (mean(trips[:, :start_lat]), mean(trips[:, :start_lon]))
)
	if !original
		percentile = min.(max.(1 .- op_cost ./ delay_cost ./ r, 0), 1)
		q = quantile.(eachrow(L_train), percentile)
		q = max.(mean(L_train, dims=2), q)
		rta_trips = copy(trips)
		rta_trips.stop_time .+= q
	
		if rta_only
			# replace trips
			trips = rta_trips
			L_train = L_train .- q
			L_test = L_test .- q
		else
			# add RTA trips
			trips = vcat(trips, rta_trips)
			L_train = vcat(L_train, L_train .- q)
			L_test = vcat(L_test, L_test .- q)
			r = vcat(r, r)
		end
	end
	L_train = vcat(zeros(Float64, 1, size(L_train, 2)), L_train)
	L_test = vcat(zeros(Float64, 1, size(L_test, 2)), L_test)

	n = size(trips, 1) + 1
	C = zeros(Float64, n, n)
    C[1, 2:end] .= veh_cost # cost per vehicle
	B = zeros(Float64, n, n)
	G = zeros(Bool, n, n)
	G[1, 2:end] .= 1 # add link from depot to each trip
	G[2:end, 1] .= 1 # add link from each trip to depot
	D = zeros(Float64, n, n) # deadhead time between start/stop points
	V = zeros(Bool, n, n)

	for i ∈ 1:n-1
		D[1, i+1] = haversine(depot_loc, (trips[i, :start_lat], trips[i, :start_lon]), 6372.8) / averageSpeed 
		D[i+1, 1] = haversine((trips[i, :stop_lat], trips[i, :stop_lon]), depot_loc, 6372.8) / averageSpeed 
		C[1, i+1] += D[1, i+1] * op_cost
		C[i+1, 1] = D[i+1, 1] * op_cost
		if !original && (i > n/2.0 || rta_only)
			C[1, i+1] += q[i-size(q, 1)*(!rta_only)] * op_cost
		end
		for j ∈ 1:n-1
            i == j && continue
			# !original && !rta_only && i % ((n-1) / 2) == j % ((n-1) / 2) && continue # duplicate trips cannot be connected
			timeDiff = trips[j, :start_time] - trips[i, :stop_time]
			timeDiff < 0 && continue
			coords_1 = (trips[i, :stop_lat], trips[i, :stop_lon])
			coords_2 = (trips[j, :start_lat], trips[j, :start_lon])
			distance = haversine(coords_1, coords_2, 6372.8)
			D[i+1, j+1] = distance / averageSpeed
			if distance / averageSpeed <= timeDiff
				G[i+1, j+1] = 1 # add link if vehicle can deadhead from i -> j
                B[i+1, j+1] = timeDiff - distance / averageSpeed
				if timeDiff < 3.0
	                C[i+1, j+1] = (distance / averageSpeed + timeDiff) * op_cost # cost per hour
				else
					V[i+1, j+1] = 1
					d1 = haversine(coords_1, depot_loc, 6372.8)
					d2 = haversine(depot_loc, coords_2, 6372.8)
					D[i+1, j+1] = (d1+d2) / averageSpeed
					C[i+1, j+1] = ((d1 + d2) / averageSpeed) * op_cost * 2 # return to depot
				end
				if !original && (j > n/2.0 || rta_only)
					C[i+1, j+1] += q[j-size(q, 1)*(!rta_only)] * op_cost
				end
			end
		end
	end

	# using longest path to tighten big-M
	g = SimpleDiGraph(G[2:end, 2:end])
	m = -1*minimum(Graphs.dijkstra_shortest_paths(g, findall(G[1,2:end]), -1*ones(Int, n-1, n-1)).dists)
	M = sum(sort(maximum.(eachrow(L_train)), rev=true)[1:m])
    B[2:end, 1] .= M # *infinite* buffer time when returning to depot

	return VSPInstance(n, M, op_cost, delay_cost, veh_cost, L_train, L_test, r, C, B, D, V, G, trips)
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