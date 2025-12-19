include("Rabi.jl")

const PATH = ""

function CalculateWigner(rabii, λf; wignerMesh=401, limits=2.0, firstIndex=1, maxt=200, numt=1000)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    println(rabii.N)

    gs = SingleWellState(rabii)	

    result = Wigner(rabif; Ψ0=gs, firstIndex=firstIndex, lastIndex=-1,
    maxt=maxt, numt=numt, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), marginals=false,
    saveData=true, saveGraph=false, showGraph=false)

    ev = ExpectationValues(rabif, [:Jx=>Jx(rabif), :Jy=>Jy(rabif), :Jz=>Jz(rabif)]; Ψ0=gs, maxt=maxt, numt=2001, asymptotics=false)
end

CalculateWigner(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -sqrt(2.0)/5, maxt=400)
