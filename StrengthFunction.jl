include("Rabi.jl")

const PATH = ""

gr()
default(size=(1920,1080))

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
δ = 0.5

n = 250

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ, j=1//2)
energies, vectors = eigenstates(rabii)
Ψ = SingleWellState(rabii)

# Level Dynamics
_, pld = LevelDynamics(rabii, ps=LinRange(-2, 2, 1001), limit=200, ylims=(-0.5,2))


# Ground state
rabif = Rabi(N=n, R=R, λ=λf, δ=δ, j=1//2)

# Strength function
sf = StrengthFunction(rabif, Ψ)

p = scatter(sf[1], sf[2], markeralpha=0.7, markerstrokewidth=0, markersize=6, xlims=(-2, 2))
display(p)

_, p = Overlap(vectors; limit=length(vectors))

p = scatter(pld, λf .- 5.0 .* sf[2], sf[1], markeralpha=0.9, markerstrokewidth=0, markersize=8)
p = plot(p, [λi, λf], [energies[1], ExpectationValue(H(rabif), Ψ) / rabif.R], linewidth=3, arrow=(:closed, 1.0), xlabel="λ", ylabel="E")
display(plot(p))
