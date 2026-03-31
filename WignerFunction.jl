include("Rabi.jl")
include("ClassicalTrajectory.jl")

const PATH = ""

gr()
default(size=(1920,1080))

function WignerFunctions(rabii; λf=0.5, wignerMesh=301, limits=1.6, maxt=30, numt=31, log=false, firstIndex=1, lastIndex=-1, kwargs...)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)	

    if log
        clim = (-6, 6)
    else
        clim = (-0.08, 0.08)
    end

    Wigner(rabif; Ψ0=gs, operators=[:Jx=>Jx(rabif), :q=>X(rabif), :Jy=>Jy(rabif), :p=>P(rabif), :Jz=>Jz(rabif), :P=>projector(gs)], operatorsColor=[:gray, :red, :gray, :red, :gray, :blue], operatorsMarker=[:diamond, :square, :diamond, :square, :diamond, :pentagon],
    operatorsLayout=(7, 2), firstIndex=firstIndex, lastIndex=lastIndex,
    maxt=maxt, numt=numt, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), 
    clim=clim, saveData=true, saveGraph=true, showGraph=true, log=log, postProcess=ClassicalToWigner, postProcessParams=rabii.λ, kwargs...)
end

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5

δ = 0.5

# n = 300
j = 1 // 2

# Initial and final systems
rabii = Rabi(R=R, λ=λi, δ=δ, j=j)

limits = 1.5
wignerMesh = 301

maxt = 30*pi
numt = 300
firstIndex = 1

WignerFunctions(rabii; λf=λf, wignerMesh=wignerMesh, limits=limits, maxt=maxt, numt=numt, log=false, firstIndex=firstIndex, lastIndex=-1)