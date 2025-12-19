include("Rabi.jl")

const PATH = ""

gr()
default(size=(1920,1080))

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
μ = 0.0

δ = 0.5

# n = 300
j = 4 // 2

# Initial and final systems
rabii = Rabi(R=R, λ=λi, δ=δ, j=j, μ=μ, ν=μ)
rabif = Copy(rabii; λ=λf)

gs = SingleWellState(rabii)	

limits = 1.5
wignerMesh = 301

clim = (-0.08, 0.08)
firstIndex = 1

Wigner(rabif; Ψ0=gs, operators=[:Jx=>Jx(rabif), :q=>X(rabif), :Jy=>Jy(rabif), :p=>P(rabif), :Jz=>Jz(rabif), :P=>projector(gs)], operatorsColor=[:gray, :red, :gray, :red, :gray, :blue], operatorsMarker=[:diamond, :square, :diamond, :square, :diamond, :pentagon],
operatorsLayout=(7, 2), firstIndex=firstIndex, lastIndex=-1,
maxt=30*pi, numt=6, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), marginals=false,
clim=clim, saveData=true, saveGraph=true, showGraph=false)