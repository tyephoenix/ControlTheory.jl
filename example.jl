using .ControlTheory
using DifferentialEquations
using Plots;gr()

# Variables
tSpan = (0, 60)

Ca_0 = 20.0
Cb_0 = 30.0
Cc_0 = 0.0
k = 0.02


# Inputs
function pump(M::Matrix{Float64})
    function pump0(p, t)
        return max(p.p,0)*M
    end
    return pump0
end
ca_pump = Pump(pump([1.0;0.0;0.0;;]))
cb_pump = Pump(pump([0.0;1.0;0.0;;]))


# True System
function reaction(du, u, p, t)
    du[1] = -k * u[1] * u[2]
    du[2] = -k * u[1] * u[2]
    du[3] = k * u[1] * u[2]
end
# USE StateModel([1.0;0.0;;], reaction, inputs=[ca_pump, cb_pump]) for a system with pumps! 
# Below represents system with no pumps attached
cstr = StateModel([1.0;0.0;0.0;;], reaction) 

# Filtering Data
# T is the delta T of the differential equation. τ is the time constant.
filter = ExponentialFilter([0.0;1.0;0.0;;], cstr, 0.01, 0.0001)

# Measurements (Ca + Cb) and (Cc)
function measurement(u)
    return [u[1] + u[2]; u[3];;]
end
sensor = MeasurementSensor(filter, measurement)

# Observer
L = [0.023 0.0; 0.53 0.026; 0.44 0.50]
# USE LuenbergerObserver([0.0;1.0;;],[sensor], reaction,L=L) to have predetermined L matrix
# Below uses optimization funcion using NLopt to find L matrix
# With a predetermined L matrix (L), there is no need for lines 52-53
luenberger = LuenbergerObserver([0.0;0.0;1.0;;],[sensor], reaction)

# Optimization
package = optimize(luenberger, tSpan, [0.0;0.0;0.0;;], 20, maxtime=10)
solve(ODEProblem(package[1], [Ca_0; Cb_0; Cc_0], tSpan), callback=package[2], tstops=package[3])

# Controller
mpc = ModelPredictiveController([ca_pump, cb_pump], luenberger, [10.0;50.0;0.0;;], [0.1 0.0 0.0; 0.0 0.1 0.0])


# System
u0 = [Ca_0 Ca_0 0.0; Cb_0 Cb_0 0.0; Cc_0 Cc_0 0.0]
problem = ODEProblem(mpc.system(), u0, tSpan, [])
sol = solve(problem, Tsit5())


# Plot
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

graph = plot(xlim=(0, Inf), title="State Controlled Bilinear Reaction", dpi=600)
plot!(graph, sol.t, sol[4,:], label="Filtered Ca")
plot!(graph, sol.t, sol[5,:], label="Filtered Cb")
plot!(graph, sol.t, sol[6,:], label="Filtered Cc")
plot!(graph, sol.t, sol[7,:], linestyle=:dash, label="Ca*")
plot!(graph, sol.t, sol[8,:], linestyle=:dash, label="Cb*")
plot!(graph, sol.t, sol[9,:], linestyle=:dash, label="Cc*")
display(graph)
Save("example")