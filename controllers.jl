abstract type StateController <: StateProcess end
abstract type OptimizableController <: StateController end

mutable struct ModelPredictiveController <: OptimizableController

    inputs::Array{StateInput}
    setPoints::Matrix{Float64}
    observer::StateObserver

    K::Matrix{Float64}
    accuracy::Float64

    system::Function

    function system(self::ModelPredictiveController)
        function system0(du, u, p, t)
            self.observer.system()(du,u,p,t)

            ex = u * self.observer.A

            mat = self.K*(self.setPoints - ex)

            for i in 1:size(self.inputs)[1]
                self.inputs[i].p = mat[i]
            end
        end
        return system0
    end

    function ModelPredictiveController(inputs, observer, setPoints, K)
        self = new(inputs, setPoints, observer, K)
        function system0()
            system(self)
        end
        self.system = system0
        return self
    end
end