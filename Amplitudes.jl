using StatsBase
using DifferentialEquations
using LinearAlgebra

include("Rabi.jl")

const PATH = ""

gr()
default(size=(1920,1080))

""" Evolution of the partial trace and related quantities """
function PartialTraceEvolution(rabii; λf=0.5, mint=0.0, maxt=50, numt=1000, log=false)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)

    spin = Array{Float64}(undef, numt, 4)
    ptt = Array{Float64}(undef, numt, 4)
    ptev = Array{Float64}(undef, numt, Int(2 * j + 1))
    ptdiag = Array{Float64}(undef, numt, Int(2 * j + 1))
    field = Array{Float64}(undef, numt, 4)

    time = @elapsed begin
        u = U(rabif, mint)
        du = U(rabif, (maxt - mint) / numt)
        Ψ = u * gs      # Initial state

        for i in 1:numt
            ρq = ptrace(Ψ, 1)
            es, vs = QuantumOptics.eigenstates(ρq)

            for j in 1:Int(2 * j + 1)
                ptev[i, j] = if log log10(real(es[j])) else real(es[j]) end
                ptdiag[i, 1] = real(ρq.data[j, j])
            end

            ptt[i, 1] = real(ρq.data[1, 2])
            ptt[i, 2] = -imag(ρq.data[1, 2])
            ptt[i, 3] = real(ρq.data[1, 1]) - 0.5
            ptt[i, 4] = sqrt(ptt[i, 1]^2 + ptt[i, 2]^2 + ptt[i, 3]^2)

            q = ExpectationValue("X", X(rabif), Ψ, rabif)
            p = ExpectationValue("P", P(rabif), Ψ, rabif)

            field[i, 1] = sqrt(8) * rabif.λ * q
            field[i, 2] = -sqrt(8) * rabif.λ * rabif.δ * p
            field[i, 3] = 1 / (2 * rabif.j) + sqrt(8) * rabif.μ * q
            field[i, 4] = sqrt(field[i, 1]^2 + field[i, 2]^2 + field[i, 3]^2)

            for j in 1:3
                field[i, j] = -field[i, j] / field[i, 4]
            end

            spin[i, 1] = ExpectationValue("Jx", Jx(rabif), Ψ, rabif)
            spin[i, 2] = ExpectationValue("Jy", Jy(rabif), Ψ, rabif)
            spin[i, 3] = ExpectationValue("Jz", Jz(rabif), Ψ, rabif)
            spin[i, 4] = sqrt(spin[i, 1]^2 + spin[i, 2]^2 + spin[i, 3]^2)

            for j in 1:3
                spin[i, j] = spin[i, j] / spin[i, 4]
            end

            Ψ = du * Ψ
        end
    end
    println(time)

    tout = LinRange(mint, maxt, numt)

    p = plot(tout, field[:, 1:3], xlabel="t", lc = [:red :blue :green], label=["bx" "by" "bz"], lw=2, alpha=0.5)
    p = plot!(p, tout, spin[:, 1:3], lc = [:red :blue :green], label=["jx" "jy" "jz"], lw=1, ls=:dash, alpha=0.5)
    display(p)
    savefig(p, "$(PATH)components_($(rabif)).png")

    if log
        p = plot(tout, ptev, ylims=(-5, 0), lc = :black, xlabel="\$t\$", label=["ρ_spin EV 1" "ρ_spin EV 2"], ls=[:dash :solid])
        q = plot(tout, ptdiag, ylims=(-5, 0), lc = :black, xlabel="\$t\$", label=["ρ_diag 1" "ρ_diag 2"], ls=[:dash :solid])
    else
        p = plot(tout, ptev, lc = :black, xlabel="\$t\$", label=["ρ_spin EV 1" "ρ_spin EV 2"], ls=[:dash :solid])
        q = plot(tout, ptdiag, lc = :black, xlabel="\$t\$", label=["ρ_diag 1" "ρ_diag 2"], ls=[:dash :solid])
    end

    display(p)
    savefig(p, "$(PATH)ptt_($(rabif)).png")
end

function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1    
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function PartialTraceEvolutionScaling(rabii; R0=10, numR=6, λi=1.5, mint=0.0, maxt=50, numt=1000)
    pt = Array{Float64}(undef, numt, 2 * numR)

    tout = LinRange(mint, maxt, numt)
    lc = []
    label = []

    colors = [:red :blue :green :purple :orange :black]

    R = R0
    for j = 1:numR
        print(R)
        print("...")

        # Initial and final systems
        rabif = Copy(rabii; λ=λf)
        gs = SingleWellState(rabii)
    
        time = @elapsed begin
            u = U(rabif, mint)
            du = U(rabif, (maxt - mint) / numt)
            Ψ = u * gs      # Initial state

            for i in 1:numt
                ρq = ptrace(Ψ, 1)
                es, vs = QuantumOptics.eigenstates(ρq)

                pt[i, 2*j - 1] = real(es[1])
                pt[i, 2*j] = real(es[2])

                Ψ = du * Ψ
            end
        end
        println(time)

        push!(lc, colors[j])
        push!(lc, colors[j])

        push!(label, "EV1 R=$R")
        push!(label, "EV2 R=$R")

        R *= 2
    end

    lc = permutedims(lc)
    label = permutedims(label)

    p = plot(tout, pt, xlabel="\$t\$", lc=lc, label=label)
    display(p)
    savefig(p, "$(PATH)ptt_scaling.png")
end

function alpha2(λi, λf)
    return ((λf - 16 * λf * λi^4 + λi * (-1 + 4 * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4)))) /
              (λi * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4)))) / 8
end

function beta2(λi, λf)
    return 0.25 * (2 + (λi + λf * (-1 + 16 * λi^4)) / (2 * λi * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4))))
end

function StrengthFunctionMagnitude(rabii; λf_min=-0.5, λf_max=1.0, λf_num=50)
    peaks = Array{Float64}(undef, λf_num, Int(2 * rabii.j + 1))
    λfs = LinRange(λf_min, λf_max, λf_num)

    for i = eachindex(λfs)
        λf = λfs[i]
        ps = StrengthFunction(rabii; λf=λf, showgraph=false, envelope_window=Int(6 * rabii.j + 2))

        for k = eachindex(ps)
            peaks[i, k] = ps[k]
        end
    end

    p = plot(λfs, peaks, xlabel="\$λ_f\$", ylabel="Strength", title="Strength function")
    display(p)

    return λfs, peaks
end

_, pld = LevelDynamics(Rabi(R=10,N=300), ps=LinRange(-1.5, 1.5, 401), limit=150, ylims=(-1,2), saveGraph=false)

rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
StrengthFunction(rabi; λf=-sqrt(2)/5, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)

""" Fast oscillations for the quasiseparated state """
function RabiOscillations()
    λf = -0.4
    λi = 1.4
    rabi = Rabi(R=20, λ=λi, δ=0.5, j=4//2)

    maxt = 200
    numt = 2
    limits = 1.5

    WignerFunctions(rabi, λf=λf, limits=limits, wignerMesh=301, maxt=maxt, numt=numt, showGraph=false, marginals=true)
end

RabiOscillations()