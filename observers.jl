abstract type StateObserver <: StateProcess end
abstract type OptimizableObserver <: StateObserver end

mutable struct LuenbergerObserver <: OptimizableObserver

    inputs::Array{StateSensor}
    A::Matrix{Float64}
    model::Function

    L::Matrix{Float64}
    accuracy::Float64

    # observe::Function
    # set::Function
    # bounds::Function
    # format::Function

    # function observe(self::LuenbergerObserver, stateEstimator::StateEstimator, ym::Vector{Float64}, y_hat::Vector{Float64}, du, u, p, t) 
    #     mat = self.L * (ym - y_hat)
    #     i = 1
    #     while (i <= size(du)[1])
    #         du[i] = du[i] + mat[i]
    #         i = i + 1
    #     end
    # end

    # function set(self::LuenbergerObserver, stateEstimator::StateEstimator, newValue::Matrix{Float64}, accuracy::Float64)
    #     self.L = newValue
    #     self.accuracy = accuracy
    # end

    # function bounds(self::LuenbergerObserver, stateEstimator::StateEstimator)
    #     testArr = stateEstimator.measurement(stateEstimator.trueProblem.u0)
    #     s = size(testArr)[1] * size(stateEstimator.trueProblem.u0)[1]
    #     return [fill(0,s),fill(1,s)]
    # end

    # function formatArr(self::LuenbergerObserver, stateEstimator::StateEstimator, x::Array)
    #     testArr = stateEstimator.measurement(stateEstimator.trueProblem.u0)
    #     formatted = reshape(x, (size(stateEstimator.trueProblem.u0)[1], size(testArr)[1]))
    #     return formatted
    # end


    system::Function

    function system(self::LuenbergerObserver)
        function system0(du,u,p,t)
            ym = nothing
            for i in 1:size(self.inputs)[1]
                sensor = self.inputs[i]  
                if (ym == nothing)
                    ym = sensor.system()(du, u, p, t)
                else 
                    ym = vcat(ym, sensor.system()(du,u,p,t))
                end
            end

            x = u * self.A
            dx = du * self.A
            B = transpose(self.A)
            yhat = nothing
            for i in 1:size(self.inputs)[1]
                sensor = self.inputs[i]  
                if (yhat == nothing)
                    yhat = sensor.measurement(x)
                else 
                    yhat = vcat(ym, sensor.measurement(x))
                end
            end

            self.model(dx, x, p, t)

            mat = self.L * (ym - yhat)
            i = 1
            while (i <= size(dx)[1])
                dx[i] = dx[i] + mat[i]
                i = i + 1
            end

            u .= x*B - (u*self.A)*B + u
            du .= dx*B - (du*self.A)*B + du
        end
        return system0
    end

    function LuenbergerObserver(A, inputs, model; L=zeros(Float64,1,1)) 
        self = new(inputs, A, model, L, 0.0)
        function system0()
            system(self)
        end
        self.system = system0
        # function observe0(stateEstimator::StateEstimator, ym::Vector{Float64}, y_hat::Vector{Float64}, du, u, p, t)
        #     observe(self, stateEstimator, ym, y_hat, du, u, p, t)
        # end
        # self.observe = observe0
        # function set0(stateEstimator::StateEstimator, newValue::Matrix{Float64}, accuracy::Float64)
        #     set(self, stateEstimator, newValue, accuracy)
        # end
        # self.set = set0
        # function bounds0(stateEstimator::StateEstimator)
        #     bounds(self, stateEstimator)
        # end
        # self.bounds = bounds0
        # function format0(stateEstimator::StateEstimator, x::Array)
        #     formatArr(self, stateEstimator, x)
        # end
        # self.format = format0
        return self
    end
end
