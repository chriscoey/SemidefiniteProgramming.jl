# wrapper to convert SDP solver into Conic solver

# To enable Conic support from an SDP solver, define, e.g.,
# ConicModel(s::CSDPSolver) = SDPtoConicBridge(SDPModel(s))

# This file is adapted from lpqp_to_conic.jl file from MathProgBase

type SDPtoConicBridge <: AbstractConicModel
    sdpmodel::AbstractSDPModel
    c
    A
    b
    constr_cones
    var_cones
end

SDPtoConicBridge(m::AbstractSDPModel) = SDPtoConicBridge(m, nothing, nothing, nothing, nothing, nothing)

export SDPtoConicBridge

numvar(m::SDPtoConicBridge) = size(m.A,2)
numconstr(m::SDPtoConicBridge) = size(m.A,1)

function getmatdim(k)
    # n*(n+1)/2 = k
    # n^2+n-2k = 0
    # (-1 + sqrt(1 + 8k))/2
    n = (-1 + sqrt(1 + 8k)) / 2
    if n * (n+1) != 2*k
        error("sdp dim not consistent")
    end
    n
end

# To transform Conic problems into SDP problems
function loadproblem!(m::SDPtoConicBridge, c, A, b, constr_cones, var_cones)
    m.c = c
    m.A = A
    m.b = b
    m.constr_cones = constr_cones
    m.var_cones = var_cones

    # Conic form        LP form
    # min  c'x          min      c'x
    #  st b-Ax ∈ K_1     st lb <= Ax <= b
    #        x ∈ K_2         l <=  x <= u

    # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated,:SDP], give up.
    bad_cones = [:ExpPrimal, :ExpDual]
    for (cone,idxs) in var_cones
        cone in bad_cones && error("Cone type $(cone) not supported")
    end

    blk = 0
    varmap = Vector{Tuple{Int,Int,Int,Float64}}[]
    for (cone,idxs) in var_cones
        # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated,:SDP], give up.
        if cone == :Free
            for i in idxs
                blk += 2
                varmap[i] = [(blk-1,1,1,1.), (blk,1,1,-1.)]
            end
        elseif cone == :Zero
            for i in idxs
                varmap[i] = []
            end
        elseif cone == :NonNeg
            for i in idxs
                blk += 1
                varmap[i] = [(blk,1,1,1.)]
            end
        elseif cone == :NonPos
            for i in idxs
                blk += 1
                varmap[i] = [(blk,1,1,-1.)]
            end
        elseif cone == :SOC
            error("not supported yet")
        elseif cone == :SOCRotated
            error("not supported yet")
        elseif cone == :SDP
            n = getmatdim(length(idxs))
            k = 0
            for i in 1:n
                for j in i:n
                    k += 1
                    blk += 1
                    varmap[idxs[k]] = [(blk,i,j,1/sqrt(2))]
                end
            end
        end
    end
    sdp = SparseSDP(maximize=false)
    constr = 0
    constrmap = Vector{Int}(length(b))
    slackmap = Vector{Int}(length(b))
    for cone,idxs in constr_cones
        if cone == :Free
            constrmap[idxs] = 0
            slackmap[idxs] = 0
        elseif cone == :SOC
            error("not supported yet")
        elseif cone == :SOCRotated
            error("not supported yet")
        else
            for idx in idxs
                constr += 1
                constrmap[idx] = constr
            end
            if cone == :Zero
                slackmap[idxs] = (0,0,0,0.)
            elseif cone == :NonNeg
                for idx in idxs
                    blk += 1
                    slackmap[idx] = (blk,1,1,-1.)
                end
            elseif cone == :NonPos
                for idx in idxs
                    blk += 1
                    slackmap[idx] = (blk,1,1,1.)
                end
            elseif cone == :SDP
                n = getmatdim(length(idxs))
                k = 0
                for i in 1:n
                    for j in i:n
                        k += 1
                        blk += 1
                        slackmap[idxs[k]] = (blk,i,j,-1.)
                    end
                end
            end
        end
    end
    for c in 1:length(b)
        if constrmap[c] != 0
            setrhs!(sdp, constr, b[c])
        end
        blk, i, j, coef = slackmap[c]
        if blk != 0
            setcon!(sdp, blk, i, j, coef)
        end
    end
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for col = 1:n
        for k in nzrange(A, col)
            c = rows[k]
            if constrmap[c] != 0 # Free constraint
                val = vals[k]
                blk, i, j, coef = varmap[col]
                setcon!(sdp, constrmap[c], i, j, val*coef)
            end
        end
    end
end

for f in [:optimize!, :status, :getsolution, :getobjval, :getvartype]
    @eval $f(model::SDPtoConicBridge) = $f(model.sdpmodel)
end

setvartype!(model::SDPtoConicBridge, vtype) = setvartype!(model.sdpmodel, vtype)

for f in methods_by_tag[:rewrap]
    @eval $f(model::SDPtoConicBridge) = $f(model.sdpmodel)
end
