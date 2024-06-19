include("nodes.jl")

abstract type StateProcess <:CT_Node end

struct StateModel <: StateProcess

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
            dx = self.A * du
            x = self.A * u
            B = transpose(self.A)

            for i in 1:size(self.inputs)[1]
                out = self.inputs.output(t)
                dx = dx + out
            end

            self.model(dx,x,p,t)

            u .= x*B - u*self.A*B + u
            du .= dx*B - du*self.A*B + du
        end
        return system0
    end

    function StateModel(A, model; inputs=[])
        self = new(A,inputs, model)
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