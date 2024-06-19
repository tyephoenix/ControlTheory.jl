module ControlTheory
    include("nodes.jl")
    # include("controllers.jl")
    include("observers.jl")
    # include("processes.jl")

    export ControlSystem

    export Pump
    export MeasurementSensor
    export StateModel

    export LuenbergerObserver

    struct ControlSystem
        initialState::Matrix{Float64}

        ControlSystem(initialState) = new(initialState)
    end
end;