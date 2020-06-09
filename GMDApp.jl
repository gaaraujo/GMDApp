# GMDApp.jl (WIP)

import Base.parse, Base.getindex
import DelimitedFiles

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

function parse(::Type{Tuple{T}}, s::AbstractString) where T
    Tuple(parse.(T, split(s, ',')))
end


io = open("test/NGA-Sub_Intraslab.csv")
nrec = countlines(io) - 4
seekstart(io)
readline(io)
n = parse(Int64, readline(io))
T = parse.(Float64, split(readline(io), ','))
readline(io)
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
close(io)


io = open("test/filter")
filterdict = Dict()
for line in eachline(io)
    k, v = split(line)
    if v != "NA"
        if k == "FaultType"
            push!(filterdict, k => parse(Int64, v))
        else
            push!(filterdict, k => parse(Tuple{Float64}, v))
        end
    end
end
close(io)

is = ones(Bool, nrec)
for (i, rec) in enumerate(records)
    for (k, v) in filterdict
        if k == "FaultType"
            rec[k] != v && (is[i] &= false)
        else
            !(first(v) ≤ rec[k] ≤ last(v)) && (is[i] &= false)
        end
    end
end


obj = DelimitedFiles.readdlm("test/example_spectra.csv", ',', Float64, skipstart=2)
Tobj = obj[:, 1]
Sobj = obj[:, 2]

