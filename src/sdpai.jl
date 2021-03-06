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

    LPblock = find(size -> (sign(size) == -1), sizeblocks)[1]
    PSDblocks = find(size -> (sign(size) == 1), sizeblocks)

    sizeblocks[LPblock] *= -1

    numcons = sizeblocks[LPblock] + sum(binomial((sizeblocks[block] + 1), 2) for block in PSDblocks)
    b = zeros(numcons)
    A = spzeros(numcons, numvars)

    bPSDs = SparseMatrixCSC{Float64,Int64}[spzeros(size, size) for size in sizeblocks]
    APSDs = SparseMatrixCSC{Float64,Int64}[spzeros(size, size) for size in sizeblocks, col in 1:numvars]

    while !eof(io)
        line = strip(readline(io))

        if !startswith(line, '*') && !startswith(line, '#') && !startswith(line, '"') && length(line) > 0
            linevec = split(line)
            col = parse(Int, linevec[1])
            block = parse(Int, linevec[2])
            i = parse(Int, linevec[3])
            j = parse(Int, linevec[4])
            v = parse(Float64, linevec[5])

            if block == LPblock
                @assert i == j
                if col == 0
                    b[i] = -v
                else
                    A[i,col] = -v
                end
            else
                if col == 0
                    bPSDs[block][i,j] = -v
                else
                    APSDs[block,col][i,j] = -v
                end
            end
        elseif startswith(line, "*INTEGER")
            break
        end
    end

    dropzeros!(A)

    function ind_svec(n::Int, i::Int, j::Int)::Int
        indmin = min(i,j)
        indmax = max(i,j)

        accum = indmax + 1 - indmin
        for k in 1:(indmin - 1)
            accum += n + 1 - k
        end
        return accum
    end

    varcones = Tuple{Symbol,Vector{Int}}[(:Free, collect(1:numvars))]

    rowend = sizeblocks[LPblock]
    concones = Tuple{Symbol,Vector{Int}}[(:NonNeg, collect(1:rowend))]

    sqrt2 = sqrt(2)

    for block in 1:numblocks
        if block == LPblock
            continue
        end

        rowendprev = rowend

        dropzeros!(bPSDs[block])
        (I,J,V) = findnz(bPSDs[block])
        for nzind in 1:length(V)
            row = rowendprev + ind_svec(sizeblocks[block], I[nzind], J[nzind])
            b[row] = (I[nzind] == J[nzind]) ? V[nzind] : V[nzind]*sqrt2
        end

        for col in 1:numvars
            dropzeros!(APSDs[block,col])
            (I,J,V) = findnz(APSDs[block,col])
            for nzind in 1:length(V)
                row = rowendprev + ind_svec(sizeblocks[block], I[nzind], J[nzind])
                A[row,col] = (I[nzind] == J[nzind]) ? V[nzind] : V[nzind]*sqrt2
            end
        end

        rowend = rowendprev + binomial((sizeblocks[block] + 1), 2)
        push!(concones, (:SDP, collect(rowendprev+1:rowend)))
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

    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, concones, varcones)
    MathProgBase.setvartype!(model, vartypes)

    return model
end
