using Gurobi
using JuMP
using LinearAlgebra
using MultiObjectiveAlgorithms
MOA = MultiObjectiveAlgorithms
using Random

gurobi_env = Gurobi.Env()
	
"""
    VSPModel

Model for the Vehicle Scheduling Problem (VSP) with propagated delays.

# Fields
- `inst::VSPInstance`: VSP instance on which the model will be applied.
- `model::JuMP.Model`: JuMP model for optimization.
- `x::Matrix{VariableRef}`: references to arc decision variables in the `model`.
- `s::Matrix{VariableRef}`: references to the propagated delay variables in the `model`.
- `endoftrip::Bool`: whether the model is optimizing over end of trip delays or not
"""
struct VSPModel
    inst::VSPInstance
    model::JuMP.Model
    x::Matrix{VariableRef}
    s::Matrix{VariableRef}
    endoftrip::Bool
end

"""
    function VSPModel(
        inst::VSPInstance[,
        L_train = nothing,
        endoftrip = true,
        warmStart = nothing,
        isInt = true,
        multiObj = false,
        silent = true,
        outputFlag = 0,
        timeLimit = 60]
    )

Create a VSP model object from `inst`.

`warmStart` can be used to initiate link variables from
an input solution.  `isInt`enforces if arc decision variables should be integer or not.
`multiObj` enforces whether the model should apply ϵ-constrained optimization on the
delay and cost objectives.
"""
function VSPModel(
    inst::VSPInstance;
    L_train = nothing,
    endoftrip = true,
    warmStart = nothing,
    isInt = true,
    multiObj = false,
    silent = true,
    outputFlag = 0,
    timeLimit = 60
)
    if multiObj
        model = Model(
            () -> MOA.Optimizer(() -> Gurobi.Optimizer(gurobi_env))
        )
        set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    else
        model = Model(() -> Gurobi.Optimizer(gurobi_env))
    end
    silent && set_silent(model)
    set_attribute(model, "OutputFlag", outputFlag)
    set_attribute(model, "TimeLimit", timeLimit)
    n = inst.n
    if isnothing(L_train)
        L_train = inst.L_train
    end
    n_train = size(L_train, 2)
    r = inst.r
    M = inst.M
    C = inst.C
    B = inst.B
    G = inst.G

    # decision variable for arc i -> j
    if isInt
        @variable(model, x[1:n, 1:n] >= 0, Bin)
    else
        @variable(model, x[1:n, 1:n] >= 0)
    end
    # variable for propagated delay at trip i
    @variable(model, d[1:n-1, 1:n_train] >= 0)
    # variable for end of trip delay at trip i
    if endoftrip
        @variable(model, s[1:n-1, 1:n_train] >= 0)
    end
    # nonlinear variable x_ij * d_i
    @variable(model, ϕ[1:n-1, 1:n-1, 1:n_train] >= 0)
    # warm start with provided solution
    if !isnothing(warmStart)
        set_start_value.(x, warmStart.x)
    end
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # variable constraints
    # nonlinear version
    # @constraint(model, [i = 1:n-1, j = 1:n_train], d[i, j] >= sum(x[2:end, i+1] .* (d[:, j] .+ L_train[2:end, j] .- B[2:end, i+1])))
    # linear version
    @constraint(model, [i = 1:n-1, j = 1:n_train], d[i, j] >= sum(ϕ[:, i, j] .+ x[2:end, i+1] .* (L_train[2:end, j] .- B[2:end, i+1])))
    # McCormick constraints for nonlinear variable
    @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train], ϕ[i, j, k] <= M * x[i+1, j+1])
    @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train], ϕ[i, j, k] <= d[i, k])
    @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train], ϕ[i, j, k] >= d[i, k] - M * (1 - x[i+1, j+1]))
    # end of trip delay
    if endoftrip
        @constraint(model, [i = 1:n-1, j = 1:n_train], s[i, j] >= d[i, j] + sum(x[:, i+1]) * L_train[i+1, j])
        # @constraint(model, [i = 1:n-1, j = 1:n_train], s[i, j] >= sum(ϕ[i, :, j]) + sum(x[:, i+1]) * L_train[i+1, j])
    end  
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed (and no duplicates)
    @constraint(model, [i = 1:Int((n-1)/2)], sum(x[:, i+1]) + sum(x[:, i+1+Int((n-1)/2)]) == 1)
    # minimize total propagated delay and link costs
    @expression(
        model,
        delay_expr,
        endoftrip ? sum(inst.delay_cost * r' * s) / n_train : sum(inst.delay_cost * r' * d) / n_train
    )
    @expression(model, cost_expr, sum(C .* x))
    if multiObj
        @objective(model, Min, [delay_expr, cost_expr])
    else
        @objective(model, Min, delay_expr + cost_expr)
    end

    return VSPModel(inst, model, x, endoftrip ? s : d, endoftrip)
end

"""
    FirstStageProblem

Model for the Bender's decomposition first-stage problem of the Vehicle Scheduling
Problem (VSP) with propagated delays.

# Fields
- `inst::VSPInstance`: VSP instance on which the model will be applied.
- `model::JuMP.Model`: JuMP model for optimization.
- `x::Matrix{VariableRef}`: references to arc decision variables in the `model`.
- `q::VariableRef`: references to the second-stage problem objective value.
"""
struct FirstStageProblem
    inst::VSPInstance # VSP instance
    model::JuMP.Model # VSP model
    x::Matrix{VariableRef} # decision variable matrix
    # q::VariableRef # second-stage objective value
    q::Vector{VariableRef} # second-stage objective values
end

"""
    FirstStageProblem(
        inst::VSPInstance[,
        L_train = nothing,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a first-stage problem object from `inst`.
"""
function FirstStageProblem(
    inst::VSPInstance;
    L_train = nothing,
    silent = true,
    outputFlag = 0,
    timeLimit = 60
)
    model = Model(
        () ->  Gurobi.Optimizer(gurobi_env)
    )
    silent && set_silent(model)
    set_attribute(model, "OutputFlag", outputFlag)
    set_attribute(model, "TimeLimit", timeLimit)
    n = inst.n
    L_train = inst.L_train
    n_train = size(L_train, 2)
    C = inst.C
    G = inst.G

    # decision variable for arc i -> j
    @variable(model, x[1:n, 1:n], Bin)
    # second-stage objective
    # @variable(model, q >= 0) # we start with a value of 0
    @variable(model, q[1:n_train] >= 0) # we start with a value of 0
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed (and no duplicates)
    @constraint(model, [i = 1:Int((n-1)/2)], sum(x[:, i+1]) + sum(x[:, i+1+Int((n-1)/2)]) == 1)
    # minimize second-stage objective and link costs
    @objective(model, Min, inst.delay_cost / n_train * sum(q) + sum(C .* x))

    return FirstStageProblem(inst, model, x, q)
end

"""
    SecondStageProblem

Model for the Bender's decomposition second-stage problem of the Vehicle Scheduling
Problem (VSP) with propagated delays.

# Fields
- `model::JuMP.Model`: JuMP model for optimization.
- `p::Vector{VariableRef}`: references to propagated delay dual variables in the `model`.
"""
struct SecondStageProblem
    model::JuMP.Model # VSP model
    p::Vector{VariableRef} # second-stage dual variable
    q::Vector{VariableRef} # second-stage dual variable
end

"""
    SecondStageProblem(
        x::Matrix{Float64},
        inst::VSPInstance,
        s::Int[,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a second-stage problem object with arc decision variables, `x`, from
the FirstStageProblem and `inst`.  `s` is the stochastic scenario.
"""
function SecondStageProblem(
    x::Matrix{Float64},
    inst::VSPInstance,
    s::Int;
    silent = true,
    outputFlag = 0,
    timeLimit = 60
)
    model = Model(
        () ->  Gurobi.Optimizer(gurobi_env)
    )
    silent && set_silent(model)
    set_attribute(model, "OutputFlag", outputFlag)
    set_attribute(model, "TimeLimit", timeLimit)
    set_attribute(model, "InfUnbdInfo", 1) # to return extreme rays
    n = inst.n-1
    L_train = inst.L_train
    B = inst.B
    r = inst.r

    # dual variable associated with each propagated delay constraint
    @variable(model, p[1:n] >= 0)
    # dual variable associated with each end of trip delay constraint
    @variable(model, q[1:n] >= 0)
    # dual constraints
    @constraint(model, [i = 1:n; sum(fs_x[i+1, :]) == 0], p[i] == 0)
    @constraint(model, [i = 1:n; sum(fs_x[i+1, :]) == 0], q[i] == 0)
    # @constraint(model, [i = 1:n], p[i] - fs_x[i, :]' * p - q[i] <= 0)
    @constraint(model, [i = 1:n], p[i] - x[i+1, 2:end]' * p - q[i] <= 0)
    @constraint(model, [i = 1:n], q[i] <= r[i])
    # dual objective
    # @objective(model, Max, sum((fs_x' * L_train .- sum((fs_x .* B)', dims=2)) .* p .+ L_train .* q .* sum(fs_x', dims=2)))
    @objective(model, Max, sum([p[i] * sum(x[2:end, i+1] .* (L_train[2:end, s] .-  B[2:end, i+1])) + q[i] * L_train[i+1, s] * sum(x[:, i+1]) for i ∈ 1:n]))

    return SecondStageProblem(model, p, q)
end

"""
    get_pq(x::Matrix{Float64}, inst::VSPInstance)

Calculate optimal dual variables and objective value for the Bender's decomposition second-stage problem
of the schedule-aware RTA model.
"""
function get_pq(x::Matrix{Float64}, inst::VSPInstance)
    L_train = inst.L_train
    r = inst.r
    B = inst.B
    x = convert(Matrix{Bool}, round.(x))
    schedules = generate_blocks(x)
    p = zeros(Float64, size(L_train))
    q = zeros(Float64, size(L_train))
    # q = L_train .> 0
    obj = 0.0

    # mask = vec(sum(x, dims=1)) .> 0
    # q = q .* mask
    # q = q .* vcat(0, r)

    for s in schedules
        t = s[1]+1
        mask = L_train[t, :] .> 0
        q[t, mask] .= r[t-1]
        obj += sum(L_train[t, :] .* q[t, :])
        length(s) == 1 && continue

        t = s[end]+1
        prev_t = s[end-1]+1
        mask1 = L_train[prev_t, :] .- B[prev_t, t] .> 0
        mask2 = L_train[prev_t, :] .- B[prev_t, t] .+ L_train[t, :] .> 0
        mask3 = L_train[t, :] .> 0
        q[t, mask1 .& mask2] .= r[t-1]
        q[t, .!mask1 .& mask3] .= r[t-1]
        p[t, mask1] .= q[t, mask1]
        obj += sum((L_train[prev_t, :] .- B[prev_t, t]) .* p[t, :]) + sum(L_train[t, :] .* q[t, :])

        for i in length(s)-1:-1:2
            t = s[i]+1
            prev_t = s[i-1]+1
            next_t = s[i+1]+1
            mask1 = L_train[prev_t, :] .- B[prev_t, t] .> 0
            mask2 = L_train[prev_t, :] .- B[prev_t, t] .+ L_train[t, :] .> 0
            mask3 = L_train[t, :] .> 0
            q[t, mask1 .& mask2] .= r[t-1]
            q[t, .!mask1 .& mask3] .= r[t-1]
            p[t, mask1] .= p[next_t, mask1] .+ q[t, mask1]
            obj += sum((L_train[prev_t, :] .- B[prev_t, t]) .* p[t, :]) + sum(L_train[t, :] .* q[t, :])
        end
    end

    # for s in schedules
    #     obj += sum(L_train[s[1]+1, :] .* q[s[1]+1, :])
    #     length(s) == 1 && continue

    #     t = s[end]+1
    #     prev_t = s[end-1]+1

    #     mask = L_train[prev_t, :] .- B[prev_t, t] .> 0
    #     p[t, mask] .= r[t-1] .+ q[t, mask]
    #     obj += sum((L_train[prev_t, :] .- B[prev_t, t]) .* p[t, :]) + sum(L_train[t, :] .* q[t, :])

    #     for i in length(s)-1:-1:2
    #         t = s[i]+1
    #         prev_t = s[i-1]+1
    #         next_t = s[i+1]+1
    #         mask = L_train[prev_t, :] .- B[prev_t, t] .> 0
    #         p[t, mask] .= p[next_t, mask] .+ q[t, mask]
    #         obj += sum((L_train[prev_t, :] .- B[prev_t, t]) .* p[t, :]) + sum(L_train[t, :] .* q[t, :])
    #     end
    # end

    return p, q, obj
end

# to measure efficiency between custom algorithm and LP
global lp_callback_runtimes = []
global get_p_callback_runtimes = []

"""
    add_benders_callback!(fs::FirstStageProblem[, silent::Bool = true])

At every integer solution of `fs`, add any violated optimality constraints from
the second-stage problem.
"""
function add_benders_callback!(fs::FirstStageProblem; silent::Bool = true, tol = 1e-4)

    # See Benders progress
	!silent && unset_silent(fs.model)

    function benders_callback(cb_data)
        
        # only run at integer solutions
		status = callback_node_status(cb_data, fs.model)
		if status != MOI.CALLBACK_NODE_STATUS_INTEGER
			return  
		end

        # generate and solve second stage problem
        S = size(fs.inst.L_train, 2)
        n = fs.inst.n
        this_x = callback_value.(cb_data, fs.x)
        Q = callback_value.(cb_data, fs.q)
        # solving via LP - can comment from here ...
        new_cuts = []
        this_obj = 0.0
        this_p = zeros(Float64, size(fs.inst.L_train))
        this_q = zeros(Float64, size(fs.inst.L_train))
        for s in 1:S
            ss = SecondStageProblem(this_x, fs.inst, s)
            optimize!(ss.model)
            p = value.(ss.p)
            this_p[2:end, s] .= p
            q = value.(ss.q)
            this_q[2:end, s] .= q
            this_obj = objective_value(ss.model)
            if this_obj > Q[s] + tol
                push!(lp_callback_runtimes, solve_time(ss.model))
                push!(new_cuts, @build_constraint(
                    sum((fs.x' * fs.inst.L_train[:, s] .- sum((fs.x .* fs.inst.B)', dims=2)) .* this_p[:, s] .+ fs.inst.L_train[:, s] .* this_q[:, s] .* sum(fs.x', dims=2)) <= fs.q[s]
                    ))
            end
        end
        # ... to here to only use algorithmic solution
        # this_p, this_q, this_obj = get_pq(this_x, inst)
        # push!(get_p_callback_runtimes, @elapsed get_pq(this_x, inst))
        # if this_obj > Q
        #     new_cut = @build_constraint(
        #         sum([sum((fs.x' * fs.inst.L_train[:, s] .- sum((fs.x .* fs.inst.B)', dims=2)) .* this_p[:, s] .+ fs.inst.L_train[:, s] .* this_q[:, s] .* sum(fs.x', dims=2)) for s ∈ 1:S]) <= fs.q
        #         )
        # else
        #     return
        # end
        length(new_cuts) > 0 && MOI.submit.(fs.model, MOI.LazyConstraint(cb_data), new_cuts)
        # MOI.submit(fs.model, MOI.LazyConstraint(cb_data), new_cut)
    end

	set_attribute(
		fs.model,
		MOI.LazyConstraintCallback(),
		benders_callback
	)

end
	
"""
    MCFModel

Model for the Vehicle Scheduling Problem (VSP).

# Fields
- `inst::VSPInstance`: VSP instance on which the model will be applied.
- `model::JuMP.Model`: JuMP model for optimization.
- `x::Matrix{VariableRef}`: references to arc decision variables in the `model`.
"""
struct MCFModel
    inst::VSPInstance # VSP instance
    model::JuMP.Model # min cost flow model
    x::Matrix{VariableRef} # decision variable matrix
end

"""
    MCFModel(
        inst::VSPInstance[,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a min-cost flow model object from `inst`.
"""
function MCFModel(
    inst::VSPInstance;
    duplicates = true,
    silent = true,
    outputFlag = 0,
    timeLimit = 60
)
    model = Model(() -> Gurobi.Optimizer(gurobi_env))
    silent && set_silent(model)
    set_attribute(model, "OutputFlag", outputFlag)
    set_attribute(model, "TimeLimit", timeLimit)
    n = inst.n
    C = inst.C
    G = inst.G

    # decision variable for arc i -> j
    @variable(model, x[1:n, 1:n] >= 0, Int)
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    if duplicates
        # no duplicates allowed
        @constraint(model, [i = 1:Int((n-1)/2)], sum(x[:, i+1]) + sum(x[:, i+1+Int((n-1)/2)]) == 1)
    else
        # require each trip is completed
        @constraint(model, [i = 2:n], sum(x[:, i]) == 1)
    end
    # minimize total link costs
    @objective(model, Min, sum(C .* x))

    return MCFModel(inst, model, x)
end

# (NOT USING)
struct DelayModel
    model::JuMP.Model # model to obtain trip delays
    s::Vector{VariableRef} # decision variable matrix
end

# (NOT USING) model to minimize s during subgradient algorithm
function DelayModel(
    λ::Vector{Float64},
    γ::Vector{Float64},
    x::Matrix{Float64},
    M::Float64;
    silent = true,
    outputFlag = 0,
    timeLimit = 60
)
    model = Model(() -> Gurobi.Optimizer(gurobi_env))
    silent && set_silent(model)
    set_attribute(model, "OutputFlag", outputFlag)
    set_attribute(model, "TimeLimit", timeLimit)
    n = size(x, 1)

    # decision variable for propagated delay
    @variable(model, 0 <= s[1:n] <= M)
    # minimize intermediate objective
    @objective(
        model,
        Min,
        sum(s .* (M .- γ .- λ) .+ x * λ)
    )

    return DelayModel(model, s)
end