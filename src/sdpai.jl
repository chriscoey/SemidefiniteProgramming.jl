using MathProgBase


function readsdpai(all_bin::Bool, io::IO)
    numvars = Int()
    numblocks = Int()
    sizeblocks = Vector{Int}()
    c = Vector{Float64}()
    indline = 0

    while !eof(io)
        line = strip(readline(io))

        if !startswith(line, '*') && !startswith(line, '#') && !startswith(line, '\"') && length(line) > 0
            indline += 1
            if indline == 1
                numvars = parse(Int, split(line)[1])
            elseif indline == 2
                numblocks = parse(Int, split(line)[1])
            elseif indline == 3
                sizeblocks = [parse(Int, t) for t in split(line)[1:numblocks]]
            elseif indline == 4
                c = [parse(Float64, t) for t in split(line)[1:numvars]]
                break
            end
        end
    end

    indLPblock = find(size -> (sign(size) == -1), sizeblocks)
    indPSDblocks = find(size -> (sign(size) == 1), sizeblocks)

    numcons = sizeblocks[indLPblock] + sum(binomial((sizeblocks[ind] + 1), 2) for ind in indPSDblocks)
    b = zeros(numcons)
    A = spzeros(numcons, numvars)

    bPSDs = SparseMatrixCSC{Float64,Int64}[spzeros(size, size) for size in sizeblocks]
    APSDs = SparseMatrixCSC{Float64,Int64}[spzeros(size, size) for size in sizeblocks, col in 1:numvars]

    while !eof(io)
        line = strip(readline(io))

        if !startswith(line, '*') && !startswith(line, '#') && !startswith(line, '"') && length(line) > 0
            linevec = split(line)
            indvar = parse(Int, linevec[1])
            indblock = parse(Int, linevec[2])
            i = parse(Int, linevec[3])
            j = parse(Int, linevec[4])
            v = parse(Float64, linevec[5])

            if indblock == indLPblock
                @assert i == j
                if indvar == 0
                    b[i] = -v
                else
                    A[i,j] = -v
                end
            else
                if indvar == 0
                    bPSDs[indblock][i,j] = -v
                else
                    APSDs[indblock,indvar][i,j] = -v
                end
            end
        elseif startswith(l, "*INTEGER")
            break
        end
    end

    var_cones = Tuple{Symbol,Vector{Int}}[(:Free, collect(1:numvars))]

    row = -sizeblocks[indLPblock]
    con_cones = Tuple{Symbol,Vector{Int}}[(:NonNeg, collect(1:row))]

    for block in 1:numblocks
        if block == indLPblock
            continue
        end

        prevrow = row
        for i in 1:sizeblocks[block], j in i:sizeblocks[block]
            row += 1

            b[row] = bPSDs[block][i,j]
            if i != j
                b[row] *= sqrt(2)
            end

            for col in 1:numvars
                A[row,col] = APSDs[block,col][i,j]
                if i != j
                    A[row,col] *= sqrt(2)
                end
            end
        end

        push!(con_cones, (:SDP, collect(prevrow+1:row)))
    end

    vartypes = fill(:Cont, numvars)

    while !eof(io)
        line = strip(readline(io))

        if startswith(line, '*')
            col = parse(Int, split(strip(line[2:end]))[1])
            vartypes[col] = all_bin ? :Bin : :Int
        end
    end

    return (c, A, b, con_cones, var_cones, var_types)
end


function loadsdpai(solver, filename::String; all_bin=false)
    open(filename, 'r') do f
        (c, A, b, con_cones, var_cones, var_types) = readsdpai(all_bin, f)
    end

    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, con_cones, var_cones)
    MathProgBase.setvartype!(model, var_types)

    return model
end
