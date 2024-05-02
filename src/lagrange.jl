
struct LagrangeRelaxation
    inst::VSPInstance # VSP instance
    λ::Array{Vector{Float64}} # Lagrange multipliers
    γ::Array{Vector{Float64}} # Lagrange multipliers
    s::Array{Vector{Float64}} # propagated trip delays
    s_adj::Array{Vector{Float64}} # adjusted propagated trip delays
    x::Array{Matrix{Float64}} # link decision values
    α::Vector{Float64} # step sizes
    β::Vector{Float64} # step sizes
    LB::Vector{Float64} # langrangian dual objective (lower bound)
    UB::Vector{Float64} # objective value with adjusted propagated trip delays (upper bound)
    opt_gap::Vector{Float64}
    best_idx::Vector{Int}
    MAX_ITER::Int
end

function LagrangeRelaxation(inst::VSPInstance)
    λ = zeros(Float64, inst.n)
    γ = zeros(Float64, inst.n)
    s = zeros(Float64, inst.n)
    # s = ones(Float64, inst.n)

    mcf_sol = solve!(MCFModel(inst))
    x = mcf_sol.x

    s_adj = feasibleDelays(x, inst.l, inst.B)

    LB = lagrangianDualObjective(
        x,
        s,
        inst.C,
        λ,
        γ,
        inst.l,
        inst.B,
        inst.M
    )
    UB = objectiveValue(
        x,
        s_adj,
        inst.C,
        inst.M
    )

    return LagrangeRelaxation(
        inst,
        [λ],
        [γ],
        [s],
        [s_adj],
        [x],
        [0.0],
        [0.0],
        [LB],
        [UB],
        [],
        [1],
        200
    )
end

function update!(
    LR::LagrangeRelaxation
)
    UB = minimum(LR.UB)
    LB = LR.LB[end]
    opt_gap = (UB - LB) / UB
    push!(LR.opt_gap, opt_gap)
    if opt_gap < 0.05
        println("Solution found within 1% optimality gap.")
        return true
    elseif size(LR.opt_gap, 1) > 10 && all(LR.opt_gap[end-10] .< LR.opt_gap[end-10:end] .+ eps())
        # println("Halted due to lack of progress.")
        # return true
    end
    n = LR.inst.n
    s = LR.s[end]
    x = LR.x[end]
    s_adj = LR.s_adj[end]
    λ = LR.λ[end]
    γ = LR.γ[end]
    l = LR.inst.l
    B = LR.inst.B
    M = LR.inst.M
    C = LR.inst.C

    α = αNew(
        UB,
        LB,
        x,
        s,
        l,
        B;
        θ = 1.0
    )
    push!(LR.α, α)

    β = βNew(
        UB,
        LB,
        s;
        θ = 1.0
    )
    push!(LR.β, β)

    λ = max.(0, λ .+ α .* [x[:, i]' * (s .+ l .- B[:, i]) - s[i] for i ∈ 1:n])
    push!(LR.λ, λ)

    if all(s .== 0)
        γ .= 0
    else
        γ = max.(0, γ .+ β .* (-s))
    end
    push!(LR.γ, γ)

    ∇s = (M .- γ .- λ) .+ x * λ
    # ∇s = (0 .- γ .- λ) .+ x * λ
    # @show ∇s
    # s -= 1 / M .* ∇s
    s -= 1 / (sum(λ .+ γ) + n) .* ∇s

    # s = max.(0, s .- 1 ./ (λ .+ γ) .* ∇s)
    s[1] = 0 # depot delay is always 0
    # delay_model = DelayModel(λ, γ, x, M)
    # s = solve!(delay_model)
    push!(LR.s, s)

    inst = VSPInstance(LR.inst, λ, s)
    model = MCFModel(inst)
    mcf_sol = solve!(model)
    x = mcf_sol.x
    push!(LR.x, x)
    
    s_adj = feasibleDelays(x, l, B)
    push!(LR.s_adj, s_adj)
    
    LB = lagrangianDualObjective(
        x,
        s,
        C,
        λ,
        γ,
        l,
        B,
        M
    )
    UB = objectiveValue(
        x,
        s_adj,
        C,
        M
    )

    if LB > maximum(LR.LB)
        LR.best_idx[1] = size(LR.s_adj, 1)
    end

    push!(LR.LB, LB)
    push!(LR.UB, UB)

    return false
end

function objectiveValue(
    x::Matrix{Float64},
    s::Vector{Float64},
    C::Matrix{Float64},
    M::Float64
)
    return M * sum(s) + sum(C .* x)
end

function lagrangianDualObjective(
    x::Matrix{Float64},
    s::Vector{Float64},
    C::Matrix{Float64},
    λ::Vector{Float64},
    γ::Vector{Float64},
    l::Vector{Float64},
    B::Matrix{Float64},
    M::Float64
)
    n = size(x, 1)
    return (
        sum(s .* (M .- γ .- λ)) + 
        sum(C .* x) + 
        sum(λ' * [x[:, i]' * (s .+ l .- B[:, i]) for i ∈ 1:n])
    )
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

function αNew(
    UB::Float64,
    LB::Float64,
    x::Matrix{Float64},
    s::Vector{Float64},
    l::Vector{Float64},
    B::Matrix{Float64};
    θ::Float64 = 1.0
)
    n = size(x, 1)
    num = θ * (UB - LB)
    denom = sum([x[:, i]' * (s .+ l .- B[:, i]) - s[i] for i ∈ 1:n].^2)

    return num / denom
end

function βNew(
    UB::Float64,
    LB::Float64,
    s::Vector{Float64};
    θ::Float64 = 1.0
)
    num = θ * (UB - LB)
    denom = sum(s.^2)

    if denom == 0
        return 0
    else
        return num / denom
    end
end