include("Rabi.jl")

const PATH = ""

function CalculatePurity(rabii, λf; maxt=200, numt=2001)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    println(rabii.N)

    Ψ0 = SingleWellState(rabii)	

    ts = LinRange(0.0, maxt, numt)
    purity = Array{Any}(undef, numt)

    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabif))
    println(time, "s")

    for (i, Ψ) in enumerate(Ψt)
        ρq = ptrace(Ψ, 2)
        purity[i] = real(tr(ρq * ρq))
    end

    p = plot(ts, purity, xlabel="\$t\$", ylabel="Purity", title="Evolution of the purity", legend=false)
    savefig(p, "$(PATH)purity_$(rabii)_$(λf).png")

    Export("$(PATH)purity_$(rabii)_$(λf)", tout, purity)
end

CalculatePurity(Rabi(R=50, λ=1.5, δ=0.5, j=1//2), 0.2, maxt=300)