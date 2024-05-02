using Gurobi
using JuMP
using LinearAlgebra
using MultiObjectiveAlgorithms
MOA = MultiObjectiveAlgorithms

gurobi_env = Gurobi.Env()
	
"""
    VSPModel

Model for the Vehicle Scheduling Problem (VSP) with propagated delays.

# Fields
- `inst::VSPInstance`: VSP instance on which the model will be applied.
- `model::JuMP.Model`: JuMP model for optimization.
- `x::Matrix{VariableRef}`: references to arc decision variables in the `model`.
- `s::Vector{VariableRef}`: references to the propagated delay variables in the `model`.
"""
struct VSPModel
    inst::VSPInstance # VSP instance
    model::JuMP.Model # VSP model
    x::Matrix{VariableRef} # decision variable matrix
    s::Vector{VariableRef} # propagated delay variable vector
end

"""
    VSPModel(
        inst::VSPInstance[,
        warmStart = false,
        isInt = false,
        multiObj = false,
        silent = true,
        outputFlag = 0,
        timeLimit = 60
        ]
    )

Create a VSP model object from `inst`.

`warmStart` will initialize model variables with the min-cost flow solution.  `isInt`
enforces if arc decision variables should be integer or not.  `multiObj` enforces
whether the model should apply ϵ-constrained optimization on the delay and cost objectives.
"""
function VSPModel(
    inst::VSPInstance;
    warmStart = false,
    isInt = false,
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
    M = inst.M
    l = inst.l
    C = inst.C
    B = inst.B
    G = inst.G
    α = sum(inst.l)

    # decision variable for arc i -> j
    if isInt
        @variable(model, x[1:n, 1:n] >= 0, Int)
    else
        @variable(model, x[1:n, 1:n] >= 0)
    end
    # variable for propagated delay at trip i
    # @variable(model, s[1:n] >= 0)
    @variable(model, s[1:n-1] >= 0)
    # nonlinear variable x_ij * s_j
    # @variable(model, phi[1:n, 1:n] >= 0)
    @variable(model, ϕ[1:n-1, 1:n-1] >= 0)
    # warm start with MCF model solution
    if warmStart
        mcf_model = MCFModel(inst)
        mcf_sol = solve!(mcf_model)
        set_start_value.(x, mcf_sol.x)
        set_start_value.(s, feasibleDelays(mcf_sol.x, l, B)[2:end])
        set_start_value.(ϕ, mcf_sol.x[2:end, 2:end] .* (feasibleDelays(mcf_sol.x, l, B)[2:end] * ones(n-1)'))
    end
    # force non-existant links to 0
    @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
    # variable constraints
    # @constraint(model, [i = 1:n], s[i] >= sum(phi[:, i] .+ x[:, i] .* (l .- B[:, i])))
    @constraint(model, [i = 1:n-1], s[i] >= sum(ϕ[:, i] .+ x[2:end, i+1] .* (l[2:end] .- B[2:end, i+1])))
    # McCormick constraints for nonlinear variable
    # @constraint(model, [i = 1:n, j = 1:n], phi[i, j] <= M * x[i, j])
    # @constraint(model, [i = 1:n, j = 1:n], phi[i, j] <= s[i])
    # @constraint(model, [i = 1:n, j = 1:n], phi[i, j] >= s[i] - M * (1 - x[i, j]))
    @constraint(model, [i = 1:n-1, j = 1:n-1], ϕ[i, j] <= M * x[i+1, j+1])
    @constraint(model, [i = 1:n-1, j = 1:n-1], ϕ[i, j] <= s[i])
    @constraint(model, [i = 1:n-1, j = 1:n-1], ϕ[i, j] >= s[i] - M * (1 - x[i+1, j+1]))
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, :]) - sum(x[:, i]) == 0)
    # require each trip is completed
    @constraint(model, [i = 2:n], sum(x[:, i]) == 1)
    # minimize total propagated delay and link costs
    @expression(model, delay_expr, sum(s))
    @expression(model, cost_expr, sum(C .* x))
    if multiObj
        @objective(model, Min, [delay_expr, cost_expr])
    else
        @objective(model, Min, α * delay_expr + cost_expr)
    end

    return VSPModel(inst, model, x, s)
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