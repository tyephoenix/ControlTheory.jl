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

function pump(M::Matrix{Float64})
    function pump0(p, t)
        return max(p.p,0)*M
    end
    return pump0
end
ca_pump = Pump(pump([1.0;0.0;0.0;;]))
cb_pump = Pump(pump([0.0;1.0;0.0;;]))


# True System
cstr = StateModel([1.0;0.0;;], reaction) 

# Observer
sensor = MeasurementSensor(cstr, measurement)
L = [0.023 0.0; 0.53 0.026; 0.44 0.50]
luenberger = LuenbergerObserver([0.0;1.0;;],[sensor],reaction,L=L)

# Controller
mpc = ModelPredictiveController([ca_pump, cb_pump], luenberger, [10.0;50.0;0.0;;], [0.0 0.0 0.0; 0.0 0.5 0.0])


# System
u0 = [Ca_0 0.0; Cb_0 0.0; Cc_0 0.0]
# problem = ODEProblem(mpc.system(), u0, tSpan, [])
# sol = solve(problem, Tsit5())

optimize(luenberger, ODEProblem(cstr.system(), u0, tSpan, []))


# Plot
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

# graph = plot(xlim=(0, Inf), ylim=(0,70), title="State Controlled Bilinear Reaction", dpi=600)
# plot!(graph, sol.t, sol[1,:], label="Ca")
# plot!(graph, sol.t, sol[2,:], label="Cb")
# plot!(graph, sol.t, sol[3,:], label="Cc")
# plot!(graph, sol.t, sol[4,:], linestyle=:dash, label="Ca*")
# plot!(graph, sol.t, sol[5,:], linestyle=:dash, label="Cb*")
# plot!(graph, sol.t, sol[6,:], linestyle=:dash, label="Cc*")

# display(graph)