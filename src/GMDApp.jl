module GMDApp

import Base: parse, getindex
using DelimitedFiles: readdlm

using Gaston

using Statistics

const NMAX = 100
const DEFAULT_INPUT = "input.txt"
const WARNING_INPUT_1 =  "
WARNING: $(DEFAULT_INPUT) found in the GMDApp home directory.
Default parameters will be used during the execution.
Modify or delete $(DEFAULT_INPUT) from the GMDApp home directory 
for user-defined input.
"
const WARNING_INPUT_2 =  "
WARNING: $(DEFAULT_INPUT) not found in the GMDApp home directory.
Introduce the name of your input file: "
const WELCOME_MESSAGE = "
GMDApp -- Ground Motion Database Search Application
Version 0.1.0 64-Bit
Copyright (c) 2020 Diego Casas, Gustavo Araujo, and Alexander Arciniegas.
MIT License.
This is a short implementation in Julia of the PEER Ground Motion Database 
Web Application Algorithm for Time Series Selection.
The app is an extension of the PEER procedure for user-defined databases. 
" 
struct Record
    RSN::Int64
    EQID::Int64
    EQName::String
    Year::Int64
    SSN::Int64
    SName::String
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
        EQName,
        Year,
        SSN,
        SName,
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
            EQName,
            Year,
            SSN,
            SName,
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


"""
    read_header(hline, nt)

Returns a dictionary with field name-column number pairs.
Currently, the output dictionary only works when the records
file has field names equal to the Record struct fields.
"""
function read_header(hline::AbstractString, nt::Int64)
    fields = Dict()
    header = split(hline, ',')
    n = length(header)
    push!(fields, "FileNameVertical" => n)
    push!(fields, "FileNameHorizontal2" => n - 1)
    push!(fields, "FileNameHorizontal1" => n - 2)
    push!(fields, "S" => (n-nt-2):(n-3))
    for i in 1:(n-nt-3)
        push!(fields, header[i] => i)
    end
    fields
end

"""
    read_records(filename)

Returns a vector of Record instances read from the given file.
Also returns the periods vector for the spectra in the records.

TO DO: Deal with trailing commas.
"""
function read_records(filename)
    open(filename) do io
        nrec = countlines(io) - 3
        seekstart(io)
        readline(io) # title
        T = parse.(Float64, split(readline(io), ','))
        fields = read_header(readline(io), length(T))
        records = Vector{Record}(undef, nrec)
        for (i, line) in enumerate(eachline(io))
            row = split(line, ',')
            records[i] = Record(
                RSN = parse(Int64, row[fields["RSN"]]),
                EQID = parse(Int64, row[fields["EQID"]]),
                EQName = row[fields["Earthquake Name"]],
                Year = parse(Int64, row[fields["Year"]]),
                SSN = parse(Int64, row[fields["SSN"]]),
                SName = row[fields["Station Name"]],
                FaultType = parse(Int64, row[fields["FaultType"]]),
                Magnitude = parse(Float64, row[fields["Magnitude"]]),
                Rjb = parse(Float64, row[fields["Rjb"]]),
                ClstD = parse(Float64, row[fields["ClstD"]]),
                Vs30 = parse(Float64, row[fields["Vs30"]]),
                PGA = parse(Float64, row[fields["PGA"]]),
                PGV = parse(Float64, row[fields["PGV"]]),
                PGD = parse(Float64, row[fields["PGD"]]),
                Dur = parse(Float64, row[fields["Dur"]]),
                S = parse.(Float64, row[fields["S"]]),
                FileNameHorizontal1 = row[fields["FileNameHorizontal1"]],
                FileNameHorizontal2 = row[fields["FileNameHorizontal2"]],
                FileNameVertical = row[fields["FileNameVertical"]],
            )
        end
        records, T
    end
end


"""
    read_filters(filename)

Returns a function that is suitable for filtering an array of Record
instances from an input filter file. The output function checks for the
Record to satisfy all the conditions given in the input filter file.
It also check for all the values of the spectrum to be positive.
"""
function read_filters(filename)
    filters = Dict()
    postfilters = Dict()
    open(filename) do io
        for line in eachline(io)
            k, v = split(line)
            if k == "Nmax"
                Nmax = v == "NA" ? NMAX : min(parse(Int64, v), NMAX)
                push!(postfilters, "Nmax" => Nmax)
                continue
            end
            if k == "ScaleFactor"
                ScaleFactor = v == "NA" ? nothing : parse(Tuple{Float64}, v)
                push!(postfilters, "ScaleFactor" => ScaleFactor)
                continue
            end
            if v != "NA"
                T = k == "FaultType" ? Int64 : Float64
                push!(filters, k => parse(Tuple{T}, v))
            end
        end
    end
    function f(rec)
        for (k, v) in filters
            # returns false at the firt false condition
            # i.e. ALL conditions must be satified
            if k == "FaultType"
                !(rec[k] in v) && return false
            else
                !(first(v) ≤ rec[k] ≤ last(v)) && return false
            end
        end
        # check that all the spectrum is positive
        !any(rec.S .≤ 0)
    end
    f, postfilters
end


"""
    read_target_spectrum(filename)

Returns the periods and spectral acceleration of the target spectrum (T, S),
and arrays for the application for weights (T_w, w).
"""
function read_target_spectrum(filename)
    open(filename) do io
        readline(io) # title
        T_w = parse.(Float64, split(readline(io), ','))
        w = parse.(Float64, split(readline(io), ','))
        data = readdlm(io, ',', Float64, skipstart = 1)
        T = data[:, 1]
        S = data[:, 2]
        T, S, T_w, w
    end
end


function loginterp(x, xp, yp)
    n = length(x)
    m = length(xp)
    y = Vector{Float64}(undef, n)
    j = 1
    first_xp = first(xp)
    first_yp = first(yp)
    last_xp = last(xp)
    last_yp = last(yp)
    for i in 1:n
        if x[i] < first_xp
            y[i] = first_yp
        elseif x[i] > last_xp
            y[i] = last_yp
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


function read_input()
    files = Dict()
    if isfile(DEFAULT_INPUT)
        filename = DEFAULT_INPUT
        print(WARNING_INPUT_1)
        println()
    else
        print(WARNING_INPUT_2)
        filename = readline()
        println()
    end
    open(filename) do io
        for line in eachline(io)
            k, v = split(line)
            push!(files, k => v)
        end
    end
    files
end


function genplot(T, S_target, T_target, S_best, S_mean, S_mean_plus_std, S_mean_minus_std)
    n = length(T)
    nrec = size(S_best)[2]
    plt = plot(T, S_best[:, 1], w = "l", legend =  "''", lw = 2, lc = :gray,
        Axes(grid = :on,
          key = "left",
          axis = "loglog",
          xrange = (T_target[1], T_target[end]),
          xlabel = "'T [s]'",
          ylabel = "'Sa [g]'",
          )
    )
    for i = 1:nrec
        plot!(T, S_best[:, i], w = "l", legend =  "''", lw = 2, lc = :gray)
    end
    plot!(T, S_target, w = "l", legend = "'Target'",  lw = 2, lc = :red)
    plot!(T, S_mean, w = "l", legend = "'Mean'",  lw = 2, lc = :black)
    plot!(T, S_mean_plus_std, w = "l", legend = "'Mean+Std'",  dt = 2, lw = 2, lc = :black)
    plot!(T, S_mean_minus_std,  w = "l", legend = "'Mean-Std'",  dt = 2, lw = 2, lc = :black)
    plt
end


function write_criteria(io, filterfile, T_w, w_orig)
    println(io, "-- Summary of GMDApp Search Criteria --")
    open(filterfile) do io2
        for line in eachline(io2)
            k, v = split(line)
            println(io, k, ":,", v)
        end
    end
    print(io, "Period Array:,")
    for v in T_w
        print(io, v, ",")
    end
    println(io)
    print(io, "Weight Array:,")
    for v in w_orig
        print(io, v, ",")
    end
    println(io)
end


function write_summary(io, records, e², f)
    println(io, "-- Summary of Metadata of Selected Records --")
    if io == stdout
        println(
            io,
            "ID, RSN, EQID, MSE, ScaleFactor, Dur, Magnitude, FaultType, Rjb, ClstD, Vs30",
        )
        for j in 1:length(records)
            rec = records[j]
            println(
                io, j, ", " ,rec.RSN, ", ", rec.EQID, ", ", round(e²[j], sigdigits = 4), ", ", 
                round(f[j], sigdigits = 2), ", ", rec.Dur, ", ",  rec.Magnitude,
                ", ", rec.FaultType, ", ", rec.Rjb, ", ", rec.ClstD, ", ", 
                rec.Vs30,
            )
        end
    else
        println(
            io,
            "ID,RSN,MSE,ScaleFactor,Dur,Earthquake Name,Year,Station Name,Magnitude,FaultType,Rjb,ClstD,Vs30,FileNameHorizontal1,FileNameHorizontal2,FileNameVertical",
        )
        for j in 1:length(records)
            rec = records[j]
            println(
                io, j, "," ,rec.RSN, ",", e²[j], ",", f[j], ",", rec.Dur, ",",
                rec.EQName, ",", rec.Year, ",", rec.SName, ",",
                rec.Magnitude, ",", rec.FaultType, ",", rec.Rjb, ",", 
                rec.ClstD, ",", rec.Vs30, ",", rec.FileNameHorizontal1, ",", 
                rec.FileNameHorizontal2, ",", rec.FileNameVertical,
            )
        end
    end
end


function write_scaled_spectra(
    io,
    records,
    T,
    S_target,
    S_mean,
    S_mean_plus_std,
    S_mean_minus_std,
    S_best,
)
    println(io, "-- Scaled Spectra used in Search & Scaling --")
    # write header
    print(io, "Period,TargetSpectrum,Mean,Mean+StdDev,Mean-StdDev")
    for rec in records
        print(io, ",RSN-$(rec.RSN)")
    end
    println(io)
    # write lines with data
    for i in 1:length(T)
        print(
            io,
            T[i],
            ",",
            S_target[i],
            ",",
            S_mean[i],
            ",",
            S_mean_plus_std[i],
            ",",
            S_mean_minus_std[i],
        )
        for j in 1:length(records)
            print(io, ",", S_best[i, j])
        end
        println(io)
    end
end


function write_unscaled_spectra(io, records, T)
    println(io, "-- Unscaled Spectra --")
    # write header
    print(io, "Period")
    for rec in records
        print(io, ",RSN-$(rec.RSN)")
    end
    println(io)
    # write lines with data
    for i in 1:length(T)
        print(io, T[i])
        for rec in records
            print(io, ",", rec.S[i])
        end
        println(io)
    end
end


function main()
    print(WELCOME_MESSAGE)
    println()
    files = read_input()

    records, T = read_records(files["flatfile"])
    n = length(T)

    criteria, postfilters = read_filters(files["filters"])
    filter!(criteria, records)
    nrec = length(records)
    if nrec == 0
        println("None of the records satisfy the given filters.")
        return
    end

    T_target, S_target, T_w, w_orig = read_target_spectrum(files["spectrum"])
    S_target = loginterp(T, T_target, S_target)
    w = loginterp(T, T_w, w_orig)
    w[(T.<first(T_target)).&(T.>last(T_target))] .= 0
    w ./= sum(w) # normalize weights

    f = Vector{Float64}(undef, nrec)
    e² = Vector{Float64}(undef, nrec)
    fmin = 0
    fmax = 9999
    ScaleFactor = postfilters["ScaleFactor"]
    if !isnothing(ScaleFactor)
        fmin = first(ScaleFactor)
        fmax = last(ScaleFactor)
    end
    for (i, rec) in enumerate(records)
        fi = exp(sum(w .* log.(S_target ./ rec.S)))
        if fi < fmin
            fi = fmin
        end
        fi < fmin && (fi = fmin)
        fi > fmax && (fi = fmax)
        f[i] = fi
        e²[i] = sum(w .* log.(S_target ./ (f[i] .* rec.S)) .^ 2)
    end

    k = postfilters["Nmax"]

    nrec = length(records)
    if nrec == 0
        println("None of the records satisfy the given Scale Factor range.")
        return
    end
    k = min(k, nrec)

    ibest = partialsortperm(e², 1:k)
    write_criteria(stdout, files["filters"], T_w, w_orig)
    write_summary(stdout, records[ibest], e²[ibest], f[ibest])

    # generate output
    S_best = Matrix{Float64}(undef, n, length(ibest))
    for (j, i) in enumerate(ibest)
        rec = records[i]
        S_best[:, j] .= f[i] .* rec.S
    end
    S_mean = Vector{Float64}(undef, n)
    S_std = Vector{Float64}(undef, n)
    for i in 1:n
        S_mean[i] = mean(S_best[i, :])
        S_std[i] = exp(std(log.(S_best[i, :])))
    end
    S_mean_plus_std = S_mean .* S_std
    S_mean_minus_std = S_mean ./ S_std
    plt = genplot(T, S_target, T_target, S_best, S_mean, S_mean_plus_std, S_mean_minus_std)

    open(files["output"], "w") do io
        write_criteria(io, files["filters"], T_w, w_orig)
        println(io)
        write_summary(io, records[ibest], e²[ibest], f[ibest])
        println(io)
        write_scaled_spectra(
            io,
            records[ibest],
            T,
            S_target,
            S_mean,
            S_mean_plus_std,
            S_mean_minus_std,
            S_best,
        )
        println(io)
        write_unscaled_spectra(io, records[ibest], T)
    end
    display(plt)
    println("Press Enter to exit...")
    readline(stdin)
end


function julia_main()
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

end # module