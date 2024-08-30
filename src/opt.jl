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
- `numScenarios::Int`: number of delay scenarios considered.
- `n_train::Int`: number of delay scenarios in the training set.
- `L_train::Matrix{Float64}`: primary delays for each trip and scenario in the training set.
- `L_test::Matrix{Float64}`: primary delays for each trip and scenario in the test set.
- `x::Matrix{VariableRef}`: references to arc decision variables in the `model`.
- `s::Matrix{VariableRef}`: references to the propagated delay variables in the `model`.
"""
struct VSPModel
    inst::VSPInstance # VSP instance
    model::JuMP.Model # VSP model
    numScenarios::Int # number of scenarios
    n_train::Int # number of training scenarios
    L_train::Matrix{Float64} # primary trip delays for each scenario in the training set
    L_test::Matrix{Float64} # primary trip delays for each scenario in the test set
    x::Matrix{VariableRef} # decision variable matrix
    s::Matrix{VariableRef} # propagated delay variable matrix
end

"""
    VSPModel(
        inst::VSPInstance[,
        numScenarios = 100,
        split = 1.0,
        randomSeed = 1,
        warmStart = nothing,
        isInt = false,
        multiObj = false,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a VSP model object from `inst`.

`numScenarios` dictates how many delay scenarios to consider.  `split` is the ratio of 
scenarios that are included in the test set.  `randomSeed` to determine
randomness of trip delays.  `warmStart` can be used to initiate link variables from
an input solution.  `isInt`enforces if arc decision variables should be integer or not.
`multiObj` enforces whether the model should apply ϵ-constrained optimization on the
delay and cost objectives.
"""
function VSPModel(
    inst::VSPInstance;
    L_train = nothing,
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
    L_test = inst.L_test
    r = inst.r
    M = inst.M
    C = inst.C
    B = inst.B
    G = inst.G

    # decision variable for arc i -> j
    if isInt
        @variable(model, x[1:n, 1:n] >= 0, Int)
    else
        @variable(model, x[1:n, 1:n] >= 0)
    end
    # variable for propagated delay at trip i
    @variable(model, d[1:n-1, 1:n_train] >= 0)
    # variable for end of trip delay at trip i
    @variable(model, s[1:n-1, 1:n_train] >= 0)
    # nonlinear variable x_ij * d_j
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
    @constraint(model, [i = 1:n-1, j = 1:n_train], s[i, j] >= d[i, j] + L_train[i+1, j])    
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed
    @constraint(model, [i = 2:n], sum(x[:, i]) == 1)
    # minimize total propagated delay and link costs
    @expression(
        model,
        delay_expr,
        sum(inst.delay_cost * r' * s) / n_train
    )
    @expression(model, cost_expr, sum(C .* x))
    if multiObj
        @objective(model, Min, [delay_expr, cost_expr])
    else
        @objective(model, Min, delay_expr + cost_expr)
    end

    return VSPModel(inst, model, size(L_train, 2) + size(L_test, 2), n_train, L_train, L_test, x, s)
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
    q::VariableRef # second-stage objective value
end

"""
    FirstStageProblem(
        inst::VSPInstance[,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a first-stage problem object from `inst`.
"""
function FirstStageProblem(
    inst::VSPInstance;
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
    M = inst.M
    C = inst.C
    G = inst.G

    # decision variable for arc i -> j
    @variable(model, x[1:n, 1:n], Bin)
    # second-stage objective
    @variable(model, q >= 0) # we start with a value of 0
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed
    @constraint(model, [i = 2:n], sum(x[:, i]) == 1)
    # minimize second-stage objective and link costs
    @objective(model, Min, M * q + sum(C .* x))

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
    p::Vector{VariableRef} # second-stage dual variables
end

"""
    SecondStageProblem(
        x::Matrix{Float64},
        inst::VSPInstance[,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a second-stage problem object with arc decision variables, `x`, from
the FirstStageProblem and `inst`.
"""
function SecondStageProblem(
    x::Matrix{Float64},
    inst::VSPInstance;
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
    n = inst.n
    l = inst.l
    B = inst.B

    # dual variable associated with each propagated delay constraint
    @variable(model, p[1:n-1] >= 0)
    # dual constraint
    @constraint(model, [i = 2:n], ((2:n .== i) .- x[i, 2:n])' * p <= 1)
    # dual objective
    @objective(model, Max, sum([p[i-1] * x[2:n, i]' * (l[2:n] .- B[2:n, i]) for i ∈ 2:n]))

    return SecondStageProblem(model, p)
end

"""
    get_p(x::Matrix{Float64}, obj::Vector{Float64}[, tol::Float64 = 1e-6])

Calculate optimal dual variables for the Bender's decomposition second-stage problem
of the Vehcile Scheduling Problem (VSP) with propagated delays.
"""
function get_p(x::Matrix{Float64}, obj::Vector{Float64}; tol::Float64 = 1e-6)
    n = size(x, 1)
    x = convert(Matrix{Bool}, round.(x))
    schedules = generate_blocks(x)
    p = zeros(n-1)

    for schedule in schedules
        this_schedule = schedule
        this_obj = obj[this_schedule]
        best = 0.0
        first = 1
        last = 1
        s = 1

        for i in 1:length(this_obj)
            curr = ((i-s+1):-1:1)' * this_obj[s:i]

            if best + tol < curr
                best = curr
                first = s
                last = i
            elseif best > tol
                p[this_schedule[first:last]] .= convert(Vector{Float64}, (last-first+1):-1:1)
                best = 0.0
                s = i + 1
            elseif curr <= tol 
                s = i + 1
            end
        end

        if best > 0
            p[this_schedule[first:last]] .= convert(Vector{Float64}, (last-first+1):-1:1)
        end
    end

    return p
end

# to measure efficiency between custom algorithm and LP
global lp_callback_runtimes = []
global get_p_callback_runtimes = []

"""
    add_benders_callback!(fs::FirstStageProblem[, silent::Bool = true])

At every integer solution of `fs`, add any violated optimality constraints from
the second-stage problem.
"""
function add_benders_callback!(fs::FirstStageProblem; silent::Bool = true)

    # See Benders progress
	!silent && unset_silent(fs.model)

    function benders_callback(cb_data)
        
        # only run at integer solutions
		status = callback_node_status(cb_data, fs.model)
		if status != MOI.CALLBACK_NODE_STATUS_INTEGER
			return  
		end

        # generate and solve second stage problem
        n = fs.inst.n
        this_x = callback_value.(cb_data, fs.x)
        this_q = callback_value(cb_data, fs.q)
        # solving via LP - can comment from here ...
        ss = SecondStageProblem(this_x, fs.inst)
        optimize!(ss.model)
        push!(lp_callback_runtimes, solve_time(ss.model))
        # ... to here to only use algorithmic solution
        this_obj = [this_x[2:n, i]' * (fs.inst.l[2:n] .- fs.inst.B[2:n, i]) for i ∈ 2:n]
        this_p = get_p(this_x, this_obj)
        push!(get_p_callback_runtimes, @elapsed get_p(this_x, this_obj))
        if sum([this_p[i-1] * this_x[2:n, i]' * (fs.inst.l[2:n] .- fs.inst.B[2:n, i]) for i ∈ 2:n]) > this_q
            new_cut = @build_constraint(sum([this_p[i-1] * fs.x[2:n, i]' * (fs.inst.l[2:n] .- fs.inst.B[2:n, i]) for i ∈ 2:n]) <= fs.q)
        else
            return
        end
        MOI.submit(fs.model, MOI.LazyConstraint(cb_data), new_cut)
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
    @variable(model, x[1:n, 1:n] >= 0)
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed
    @constraint(model, [i = 2:n], sum(x[:, i]) == 1)
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