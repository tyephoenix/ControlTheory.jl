using .ControlTheory
using DifferentialEquations
using Plots;gr()

# Variables
tSpan = (0, 60)

Ca_0 = 20.0
Cb_0 = 30.0
Cc_0 = 0.0
k = 0.02


# Tank
function reaction(du, u, p, t)
    du[1] = -k * u[1] * u[2]
    du[2] = -k * u[1] * u[2]
    du[3] = k * u[1] * u[2]
end

function measurement(u)
    return [u[1] + u[2]; u[3];;]
end

u0 = [Ca_0 0.0; Cb_0 0.0; Cc_0 0.0]

cstr = StateModel([1.0;0.0;;], reaction)
sensor = MeasurementSensor(cstr, measurement)
L_nlopt = [0.023 0.0; 0.53 0.026; 0.44 0.50]
luenberger = LuenbergerObserver([0.0;1.0;;],[sensor],reaction,L=L_nlopt)

problem = ODEProblem(luenberger.system(), u0, tSpan, [])
sol = solve(problem, Tsit5())

plot(sol)