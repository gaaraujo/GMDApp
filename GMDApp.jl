# GMDApp.jl (WIP)

import Base: parse, getindex
using DelimitedFiles: readdlm


struct Record
    RSN::Int64
    EQID::Int64
    SSN::Int64
    FaultType::Int64
    Magnitude::Float64
    Rjb::Float64
    ClstD::Float64
    Vs30::Float64
    PGA::Float64
    PGV::Float64
    PGD::Float64
    Dur::Float64
    S::Vector{Float64}
    FileNameHorizontal1::String
    FileNameHorizontal2::String
    FileNameVertical::String
    function Record(;
        RSN,
        EQID,
        SSN,
        FaultType,
        Magnitude,
        Rjb,
        ClstD,
        Vs30,
        PGA,
        PGV,
        PGD,
        Dur,
        S,
        FileNameHorizontal1,
        FileNameHorizontal2,
        FileNameVertical,
    )
        new(
            RSN,
            EQID,
            SSN,
            FaultType,
            Magnitude,
            Rjb,
            ClstD,
            Vs30,
            PGA,
            PGV,
            PGD,
            Dur,
            S,
            FileNameHorizontal1,
            FileNameHorizontal2,
            FileNameVertical,
        )
    end
end


function getindex(rec::Record, k::AbstractString)
    getproperty(rec, Symbol(k))
end


function parse(::Type{Tuple{T}}, s::AbstractString) where {T}
    Tuple(parse.(T, split(s, ',')))
end


function read_records(filename)
    open(filename) do io
        nrec = countlines(io) - 3
        seekstart(io)
        readline(io) # title
        T = parse.(Float64, split(readline(io), ','))
        readline(io) # header
        records = Array{Record,1}(undef, nrec)
        for (i, line) in enumerate(eachline(io))
            row = split(line, ',')
            records[i] = Record(
                RSN = parse(Int64, row[1]),
                EQID = parse(Int64, row[2]),
                SSN = parse(Int64, row[3]),
                FaultType = 999,
                Magnitude = parse(Float64, row[4]),
                Rjb = parse(Float64, row[8]),
                ClstD = parse(Float64, row[7]),
                Vs30 = parse(Float64, row[13]),
                PGA = parse(Float64, row[20]),
                PGV = parse(Float64, row[21]),
                PGD = parse(Float64, row[22]),
                Dur = parse(Float64, row[17]),
                S = parse.(Float64, row[23:133]),
                FileNameHorizontal1 = row[134],
                FileNameHorizontal2 = row[135],
                FileNameVertical = row[136],
            )
        end
        records, T
    end
end


function read_filters(filename)
    filters = Dict()
    open("test/filter") do io
        for line in eachline(io)
            k, v = split(line)
            if v != "NA"
                if k == "FaultType"
                    push!(filters, k => parse(Int64, v))
                else
                    push!(filters, k => parse(Tuple{Float64}, v))
                end
            end
        end
    end
    function f(rec)
        for (k, v) in filters
            # returns false at the firt false condition
            # i.e. ALL conditions must be satified
            if k == "FaultType"
                rec[k] != v && return false
            else
                !(first(v) ≤ rec[k] ≤ last(v)) && return false
            end
        end
        true
    end
end


function read_target_spectrum(filename)
    data = readdlm(filename, ',', Float64, skipstart = 2)
    T = data[:, 1]
    S = data[:, 2]
    T, S
end


function loginterp(x, xp, yp)
    n = length(x)
    m = length(xp)
    y = Vector{Float64}(undef, n)
    j = 1
    for i in 1:n
        if x[i] < xp[begin]
            y[i] = yp[begin]
        elseif x[i] > xp[end]
            y[i] = yp[end]
        else
            while j < m
                if xp[j] ≤ x[i] ≤ xp[j+1]
                    k = log(yp[j+1] / yp[j]) / log(xp[j+1] / xp[j])
                    y[i] = exp(k * log(x[i] / xp[j]) + log(yp[j]))
                    # The next point could be in the same bracket,
                    # so we break the loop without accumulating j.
                    break
                end
                j += 1
            end
        end
    end
    y
end

# function main()

records, T = read_records("test/NGA-Sub_Intraslab.csv")

n = length(T)

criteria = read_filters("test/filter")
filter!(criteria, records)

nrec = length(records)

T_target, S_target = read_target_spectrum("test/example_spectra.csv")
S_target = loginterp(T, T_target, S_target)

T_w = [0.01, 1, 10, 20]
w = [1, 1, 1, 1]
w = loginterp(T, T_w, w)
w ./= sum(w)

f = Vector{Float64}(undef, nrec)
e² = Vector{Float64}(undef, nrec)
for (i, rec) in enumerate(records)
    f[i] = exp(sum(w .* log.(S_target ./ rec.S)))
    e²[i] = sum(w .* log.(S_target ./ (f[i] .* rec.S)).^2)
end

using Plots

k = 20
ibest = partialsortperm(e², 1:k)
plt = plot(T, S_target, xaxis = :log, yaxis = :log, legend = false)
for i in ibest
    rec = records[i]
    println("RSN = $(rec.RSN), e² = $(e²[i])")
    plot!(plt, T, f[i] * rec.S, color = "gray")
end
plot!(plt, T, S_target, color = "red", linewidth = 5)
# end