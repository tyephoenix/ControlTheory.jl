abstract type CT_Node end


abstract type StateInput <: CT_Node end
mutable struct Pump <: StateInput
    output::Function
    Pump(output) = new(output)
end
 

abstract type StateProcess <:CT_Node end
mutable struct StateModel <: StateProcess

    inputs::Vector{StateInput}
    A::Matrix{Float64}
    model::Function

    state::Function
    system::Function

    function state(self::StateModel, u0, t)
        system = self.system(self)

        problem = ODEProblem(system, u0, (0,t))
        solution = solve(problem, Tsit5(), saveat=t)
        return solution(t)
    end

    function system(self::StateModel)
        function system0(du,u,p,t)
            dx = du * self.A
            x = u * self.A
            B = transpose(self.A)

            for i in 1:size(self.inputs)[1]
                out = self.inputs.output(t)
                dx = dx + out
            end

            self.model(dx,x,p,t)

            u .= x*B - (u*self.A)*B + u
            du .= dx*B - (du*self.A)*B + du
        end
        return system0
    end

    function StateModel(A, model; inputs=StateInput[])
        self = new(inputs, A, model)
        function state0(u0, t)
            state(self, u0, t)
        end
        self.state = state0
        function system0()
            system(self)
        end
        self.system = system0
        return self
    end
end

abstract type StateSensor <: CT_Node end
mutable struct MeasurementSensor <: StateSensor

    parent::StateProcess
    measurement::Function

    measure::Function
    system::Function

    function measure(self::MeasurementSensor, t::Float64)
        return self.measurement(parent.state(t))
    end

    function system(self::MeasurementSensor)
        function system0(du,u,p,t)
            self.parent.system()(du,u,p,t)
            return self.measurement(u)
        end
        return system0
    end

    function MeasurementSensor(parent, measurement) 
        self = new(parent, measurement)
        function measure0(t::Float64)
            measure(self, t)
        end
        self.measure = measure0
        function system0()
            system(self)
        end
        self.system = system0
        return self
    end
end