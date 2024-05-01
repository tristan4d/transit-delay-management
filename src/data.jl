using DataFrames
using Distances
using Graphs
using Random

struct VSPInstance
	n::Int # number of trips + depot
	M::Float64 # big-M constraint variable for propagated delay
	l::Vector{Float64} # expected delay for all trips
	C::Matrix{Float64} # adjacency-cost lists for all trips
	B::Matrix{Float64} # buffer time between all trips
	G::Matrix{Bool} # connections between trips (G[i,j] = 1 if arc i -> j in G)
	trips::DataFrame # trips dataframe
end

function VSPInstance(
	trips::DataFrame;
    l::Union{Vector{Float64}, Nothing} = nothing,
	averageSpeed::Float64 = 30.0,
	randomSeed = 1,
)
	Random.seed!(randomSeed)
	n = size(trips, 1) + 1
	if isnothing(l)
        l = zeros(Float64, n)
        l[1] = 0
        # generate random delays from 0-100% of trip length
        l[2:end] .= rand(0:100, n-1) ./ 100 .* (trips[:, :stop_time] .- trips[:, :start_time])
        # l[2:end] .= rand(0:30, n-1) ./ 100 .* (trips[:, :stop_time] .- trips[:, :start_time])
    end
    M = sum(l)
	C = zeros(Float64, n, n)
    C[1, 2:end] .= 100 # cost per vehicle
	B = zeros(Float64, n, n)
	G = zeros(Bool, n, n)
    B[2:end, 1] .= M # *infinite* buffer time when returning to depot
	G[1, 2:end] .= 1 # add link from depot to each trip
	G[2:end, 1] .= 1 # add link from each trip to depot

	for i ∈ 1:n-1
		for j ∈ 1:n-1
            i == j && continue
			timeDiff = trips[j, :start_time] - trips[i, :stop_time]
			timeDiff < 0 && continue
			coords_1 = (trips[i, :stop_lat], trips[i, :stop_lon])
			coords_2 = (trips[j, :start_lat], trips[j, :start_lon])
			distance = haversine(coords_1, coords_2, 6372.8)
			if distance / averageSpeed <= timeDiff
				G[i+1, j+1] = 1 # add link if vehicle can deadhead from i -> j
				# ensure rational data
                B[i+1, j+1] = round(timeDiff - distance / averageSpeed; digits = 2)
                C[i+1, j+1] = round(distance + timeDiff; digits = 2)
			end
		end
	end

	# using longest path to tighten big-M
	# g = SimpleDiGraph(G[2:end, 2:end])
	# m = -1*minimum(Graphs.dijkstra_shortest_paths(g, findall(G[1,2:end]), -1*ones(Int, n-1, n-1)).dists)
	# M = sum(sort(l, rev=true)[1:m+1])
    # B[2:end, 1] .= M # *infinite* buffer time when returning to depot

	return VSPInstance(n, M, l, C, B, G, trips)
end

# constructor for subgradient method to generate synthetic network with adjusted link costs
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