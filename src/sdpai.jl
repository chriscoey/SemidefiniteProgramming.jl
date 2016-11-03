# Reads a SDPA file to create a MathProgBase conic model
# Detects an LP block but does not currently detect SOC or SOCRotated cones
# SDPA sparse format modified for integer constraints described at
#   http://www.opt.tu-darmstadt.de/scipsdp/downloads/data_format.txt

using MathProgBase
using GZip

function readsdpai(all_bin::Bool, io::IO)
    numvars = 0
    numblocks = 0
    sizeblocks = Vector{Int}()
    c = Vector{Float64}()
    indline = 0

    while !eof(io)
        line = strip(readline(io))

        if !startswith(line, '*') && !startswith(line, '#') && !startswith(line, '\"') && length(line) > 0
            indline += 1
            spline = split(line)
            if indline == 1
                numvars = parse(Int, spline[1])
            elseif indline == 2
                numblocks = parse(Int, spline[1])
            elseif indline == 3
                sizeblocks = [parse(Int, spline[block]) for block in 1:numblocks]
            elseif indline == 4
                c = [parse(Float64, spline[col]) for col in 1:numvars]
                break
            end
        end
    end

    indLPblock = find(size -> (sign(size) == -1), sizeblocks)[1]
    indPSDblocks = find(size -> (sign(size) == 1), sizeblocks)

    sizeblocks[indLPblock] *= -1

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
                    A[i,indvar] = -v
                end
            else
                if indvar == 0
                    bPSDs[indblock][i,j] = -v
                else
                    APSDs[indblock,indvar][i,j] = -v
                end
            end
        elseif startswith(line, "*INTEGER")
            break
        end
    end

    varcones = Tuple{Symbol,Vector{Int}}[(:Free, collect(1:numvars))]

    row = sizeblocks[indLPblock]
    concones = Tuple{Symbol,Vector{Int}}[(:NonNeg, collect(1:row))]

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

        push!(concones, (:SDP, collect(prevrow+1:row)))
    end

    vartypes = fill(:Cont, numvars)

    while !eof(io)
        line = strip(readline(io))

        if startswith(line, '*')
            col = parse(Int, split(strip(line[2:end]))[1])
            vartypes[col] = all_bin ? :Bin : :Int
        end
    end

    return (c, A, b, concones, varcones, vartypes)
end


function loadsdpai(solver, filename::String; all_bin=false)
    if endswith(filename, ".gz")
        fd = gzopen(filename, "r")
    else
        fd = open(filename, "r")
    end
    (c, A, b, concones, varcones, vartypes) = readsdpai(all_bin, fd)
    close(fd)

    @show c
    @show A
    @show b
    @show concones
    @show varcones
    @show vartypes

    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, concones, varcones)
    MathProgBase.setvartype!(model, vartypes)

    return model
end
