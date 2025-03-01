using Distributed
Distributed.@everywhere begin
    using Gurobi
    using JuMP
    using LinearAlgebra
    using MultiObjectiveAlgorithms
    MOA = MultiObjectiveAlgorithms
    using Random

    gurobi_env = Gurobi.Env()
end

Distributed.@everywhere begin
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
        x::JuMP.Containers.SparseAxisArray{JuMP.VariableRef, 2, Tuple{Int, Int}}
        s::Matrix{VariableRef}
        endoftrip::Bool
    end
end

Distributed.@everywhere begin
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
        max_vehicles = nothing,
        orig_thresh = nothing,
        warmStart = nothing,
        isInt = true,
        multiObj = false,
        silent = true,
        outputFlag = 1,
        timeLimit = 60,
        logFile = nothing,
        method = -1
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
        set_optimizer_attribute(model, "Method", method)
        if !isnothing(logFile)
            set_optimizer_attribute(model, "LogFile", logFile)
        end
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
        rta_only = inst.rta_only
        rta_mask = inst.rta_mask

        # decision variable for arc i -> j
        if isInt
            @variable(model, x[i = 1:n, j = 1:n; G[i, j]] >= 0, Bin)
        else
            @variable(model, x[i = 1:n, j = 1:n; G[i, j]] >= 0)
        end
        # variable for propagated delay at trip i
        @variable(model, d[1:n-1, 1:n_train] >= 0)
        # variable for end of trip delay at trip i
        if endoftrip
            @variable(model, s[1:n-1, 1:n_train] >= 0)
        end
        # nonlinear variable x_ij * d_i
        @variable(model, ϕ[i = 1:n-1, j = 1:n-1, k in 1:n_train; G[i+1, j+1]] >= 0)
        # warm start with provided solution
        if !isnothing(warmStart)
            for (i, j) in Tuple.(findall(G))
                set_start_value(x[i, j], warmStart.x[i, j])
            end
        end
        # force non-existant links to 0
        # @constraint(model, [i = 1:n, j = 1:n; !G[i, j]], x[i, j] == 0)
        # @constraint(model, [i = 1:n-1, j = 1:n-1; !G[i+1, j+1]], ϕ[i, j, :] == 0)
        # variable constraints
        # nonlinear version
        # @constraint(model, [i = 1:n-1, j = 1:n_train], d[i, j] >= sum(x[2:end, i+1] .* (d[:, j] .+ L_train[2:end, j] .- B[2:end, i+1])))
        # linear version
        @constraint(model, [j = 1:n-1, k = 1:n_train], d[j, k] >= sum(ϕ[i, j, k] + x[i+1, j+1] * (L_train[i+1, k] - B[i+1, j+1]) for i in 1:n-1 if G[i+1, j+1]))
        # McCormick constraints for nonlinear variable
        @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train; G[i+1, j+1]], ϕ[i, j, k] <= M * x[i+1, j+1])
        @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train; G[i+1, j+1]], ϕ[i, j, k] <= d[i, k])
        @constraint(model, [i = 1:n-1, j = 1:n-1, k = 1:n_train; G[i+1, j+1]], ϕ[i, j, k] >= d[i, k] - M * (1 - x[i+1, j+1]))
        # end of trip delay
        if endoftrip
            @constraint(model, [j = 1:n-1, k = 1:n_train], s[j, k] >= d[j, k] + sum(x[i, j+1] for i in 1:n if G[i, j+1]) * L_train[j+1, k])
            # @constraint(model, [i = 1:n-1, j = 1:n_train], s[i, j] >= sum(ϕ[i, :, j]) + sum(x[:, i+1]) * L_train[i+1, j])
        end  
        # flow constraint
        @constraint(model, [i = 1:n], sum(x[i, j] for j in 1:n if G[i, j]) - sum(x[j, i] for j in 1:n if G[j, i]) == 0)
        # require each trip is completed (and no duplicates)
        if !isnothing(max_vehicles)
            @constraint(model, sum(x[1, j] for j in 1:n if G[1, j]) <= max_vehicles)
        end
        if rta_only
            @constraint(model, [i = 1:n-1], sum(x[j, i+1] for j in 1:n if G[j, i+1]) == 1)
        else
            if !isnothing(orig_thresh)
                @constraint(model, sum(x[i+1, j+1] for i in 1:n-1-sum(rta_mask) for j in 1:n-1-sum(rta_mask) if G[i+1, j+1]) >= orig_thresh)
            end
            @constraint(model, [i = 1:n-1-sum(rta_mask); !rta_mask[i]], sum(x[j, i+1] for j in 1:n if G[j, i+1]) == 1)
            for (i, rta) in enumerate(rta_mask)
                if rta
                    paired_trip = n-sum(rta_mask[i+1:end])
                    @constraint(model, sum(x[j, i+1] for j in 1:n if G[j, i+1]) + sum(x[j, paired_trip] for j in 1:n if G[j, paired_trip]) == 1)
                end
            end
        end
        # minimize total propagated delay and link costs
        @expression(
            model,
            delay_expr,
            endoftrip ? sum(inst.delay_cost * r' * s) / n_train : sum(inst.delay_cost * r' * d) / n_train
        )
        @expression(model, cost_expr, sum(C[i, j] * x[i, j] for (i, j) in Tuple.(findall(G))))
        if multiObj
            @objective(model, Min, [delay_expr, cost_expr])
        else
            @objective(model, Min, delay_expr + cost_expr)
        end

        return VSPModel(inst, model, x, endoftrip ? s : d, endoftrip)
    end
end

Distributed.@everywhere begin
    function get_current_solution_callback!(mod::VSPModel, filename, folder)
        inst = mod.inst
        function current_solution_callback(cb_data, cb_where::Cint)
            status = callback_node_status(cb_data, mod.model)
            if status != MOI.CALLBACK_NODE_STATUS_INTEGER
                return  
            end

            runtime = Ref{Cdouble}()
            obj = Ref{Cdouble}()
            bound = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtime)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST, obj)
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBND, bound)

            gap = abs((obj[] - bound[]) / obj[]) * 100

            if gap > 5
                return
            end

            Gurobi.load_callback_variable_primal(cb_data, cb_where)

            this_x = zeros(Bool, inst.n, inst.n)
            for index in eachindex(mod.x)
                this_x[index...] = round(callback_value(cb_data, mod.x[index]))
            end

            jldsave(
                joinpath(folder, filename * "_tmp.jld2"),
                inst=inst,
                x=this_x,
                sol_time=runtime,
                objective_value=obj,
                status=status,
                gap=gap
            )
        end
        MOI.set(
            mod.model,
            Gurobi.CallbackFunction(),
            current_solution_callback
        )
    end
end

Distributed.@everywhere begin
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
end

Distributed.@everywhere begin
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
        warmStart = nothing,
        force = nothing,
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
        # warm start with provided solution
        if !isnothing(warmStart)
            set_start_value.(x, warmStart)
        end
        if !isnothing(force)
            @constraint(model, x == force)
        end
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
end

Distributed.@everywhere begin
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
end

Distributed.@everywhere begin
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
        timeLimit = 60,
        po = nothing
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
        @constraint(model, [i = 1:n; sum(x[i+1, :]) == 0], p[i] == 0)
        @constraint(model, [i = 1:n; sum(x[i+1, :]) == 0], q[i] == 0)
        # @constraint(model, [i = 1:n], p[i] - fs_x[i, :]' * p - q[i] <= 0)
        @constraint(model, [i = 1:n], p[i] - x[i+1, 2:end]' * p - q[i] <= 0)
        @constraint(model, [i = 1:n], q[i] <= r[i])
        # dual objective
        if isnothing(po)
            @objective(model, Max, sum([p[i] * sum(x[2:end, i+1] .* (L_train[2:end, s] .-  B[2:end, i+1])) + q[i] * L_train[i+1, s] * sum(x[:, i+1]) for i ∈ 1:n]))
        else
            @constraint(model, sum([p[i] * sum(x[2:end, i+1] .* (L_train[2:end, s] .-  B[2:end, i+1])) + q[i] * L_train[i+1, s] * sum(x[:, i+1]) for i ∈ 1:n]) == po)
            xc = ones(Float32, n+1, n+1) ./ 2
            @objective(model, Max, sum([p[i] * sum(xc[2:end, i+1] .* (L_train[2:end, s] .-  B[2:end, i+1])) + q[i] * L_train[i+1, s] * sum(xc[:, i+1]) for i ∈ 1:n]))
        end
        
        return SecondStageProblem(model, p, q)
    end
end

Distributed.@everywhere begin
    function solve_second_stage_problem(this_x, fs, s, Q; tol = 1e-5)
        ss = SecondStageProblem(this_x, fs.inst, s)
        optimize!(ss.model)
        this_obj = objective_value(ss.model)
        if this_obj > Q + tol
            sspo = SecondStageProblem(this_x, fs.inst, s; po = this_obj)
            optimize!(sspo.model)
            p = value.(sspo.p)
            this_p = vcat(0, p)
            q = value.(sspo.q)
            this_q = vcat(0, q)
            return @build_constraint(
                sum((fs.x' * fs.inst.L_train[:, s] .- sum((fs.x .* fs.inst.B)', dims=2)) .* this_p .+ fs.inst.L_train[:, s] .* this_q .* sum(fs.x', dims=2)) <= fs.q[s]
            )
        else
            return nothing
        end
    end
end

# to measure efficiency between custom algorithm and LP
global lp_callback_runtimes = []
global lp_num_cuts = 0

"""
    add_benders_callback!(fs::FirstStageProblem[, silent::Bool = true])

At every integer solution of `fs`, add any violated optimality constraints from
the second-stage problem.
"""
function add_benders_callback!(fs::FirstStageProblem; parallel = true, silent::Bool = true, tol = 1e-5)

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
        this_x = callback_value.(cb_data, fs.x)
        Q = callback_value.(cb_data, fs.q)
        if parallel
            args = [(s, Q[s]) for s in 1:S]
            new_cuts = Distributed.pmap(arg -> solve_second_stage_problem(this_x, fs, arg[1], arg[2]; tol = tol), args)
            new_cuts = filter(!isnothing, new_cuts)
        else
            new_cuts = []
            for s in 1:S
                ss = SecondStageProblem(this_x, fs.inst, s)
                optimize!(ss.model)
                push!(lp_callback_runtimes, solve_time(ss.model))
                this_obj = objective_value(ss.model)
                if this_obj > Q[s] + tol
                    sspo = SecondStageProblem(this_x, fs.inst, s; po = this_obj)
                    optimize!(sspo.model)
                    p = value.(sspo.p)
                    this_p = vcat(0, p)
                    q = value.(sspo.q)
                    this_q = vcat(0, q)
                    global lp_num_cuts += 1
                    push!(new_cuts, @build_constraint(
                        sum((fs.x' * fs.inst.L_train[:, s] .- sum((fs.x .* fs.inst.B)', dims=2)) .* this_p .+ fs.inst.L_train[:, s] .* this_q .* sum(fs.x', dims=2)) <= fs.q[s]
                        ))
                end
            end
        end

        length(new_cuts) > 0 && MOI.submit.(fs.model, MOI.LazyConstraint(cb_data), new_cuts)
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
    x::JuMP.Containers.SparseAxisArray{JuMP.VariableRef, 2, Tuple{Int, Int}} # decision variable matrix
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
    orig_thresh = nothing,
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
    rta_mask = inst.rta_mask

    # decision variable for arc i -> j
    @variable(model, x[i = 1:n, j = 1:n; G[i, j]] >= 0, Bin)
    # force non-existant links to 0
    # @constraint(model, [i = 1:n, j = 1:n; (i, j) in G], x[i, j] == 0)
    # flow constraint
    @constraint(model, [i = 1:n], sum(x[i, j] for j in 1:n if G[i, j]) - sum(x[j, i] for j in 1:n if G[j, i]) == 0)
    if duplicates
        if !isnothing(orig_thresh)
            @constraint(model, sum(x[i+1, j+1] for i in 1:n-1-sum(rta_mask) for j in 1:n-1-sum(rta_mask) if G[i+1, j+1]) >= orig_thresh)
        end
        # no duplicates allowed
        @constraint(model, [i = 1:n-1-sum(rta_mask); !rta_mask[i]], sum(x[j, i+1] for j in 1:n if G[j, i+1]) == 1)
        for (i, rta) in enumerate(rta_mask)
            if rta
                paired_trip = n-sum(rta_mask[i+1:end])
                @constraint(model, sum(x[j, i+1] for j in 1:n if G[j, i+1]) + sum(x[j, paired_trip] for j in 1:n if G[j, paired_trip]) == 1)
            end
        end
        # for (i, rta) in enumerate(rta_mask)
        #     if rta
        #         paired_trip = n-sum(rta_mask[i+1:end])
        #         @constraint(model, sum(x[j, i+1] for j in 1:n if G[j, i+1]) + sum(x[j, paired_trip] for j in 1:n if G[j, paired_trip]) == 1)
        #     else
        #         @constraint(model, sum(x[j, i+1] for j in 1:n if G[j, i+1]) == 1)
        #     end
        # end
    else
        # require each trip is completed
        @constraint(model, [i = 2:n], sum(x[j, i] for j in 1:n if G[j, i]) == 1)
    end
    # minimize total link costs
    @objective(model, Min, sum(C[i, j] * x[i, j] for (i, j) in Tuple.(findall(G))))

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