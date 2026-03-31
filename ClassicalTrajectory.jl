using DifferentialEquations

include("Rabi.jl")

function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1   
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function Trajectory(rabif, λi, maxt, numt)
    # Initial state
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]

    timeInterval = (0.0, maxt)
    solver = TsitPap8()
    tolerance = 1E-8
    fnc = ODEFunction(EquationOfMotion!)

    result = []

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(0, maxt, step=maxt / numt))
        solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        if length(result) == 0
            push!(result, solution.t)
        end

        push!(result, [v[1] for v in solution.u])
        push!(result, [v[2] for v in solution.u])
    end

    return result
end

function ClassicalToWigner(rabif, t, p, λi)
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]
    
    timeInterval = (0.0, t)
    solver = TsitPap8()
    tolerance = 1E-8
    fnc = ODEFunction(EquationOfMotion!)

    colours = palette(:auto)
    colourIndex = 3

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(max(t - 10, 0), t, step=0.05))
        time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        println("m = $m, time = $(time)s, solution = $(length(solution))")

        p = plot!(p, solution, idxs=(1,2), lw=2, alpha=0.4, c=colours[colourIndex], label="m=$m", xlabel=raw"$q$", ylabel=raw"$p$")
        p = scatter!(p, [solution(t)[1]], [solution(t)[2]], markersize=10, c=colours[colourIndex], label=nothing)

        colourIndex += 1
    end

    return p
end
