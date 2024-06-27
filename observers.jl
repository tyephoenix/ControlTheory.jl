abstract type StateObserver <: StateProcess end
abstract type OptimizableObserver <: StateObserver end

mutable struct LuenbergerObserver <: OptimizableObserver

    inputs::Array{StateSensor}
    A::Matrix{Float64}
    model::Function

    L::Matrix{Float64}
    accuracy::Float64

    set::Function
    system::Function

    function set(self::OptimizableObserver, M::Matrix{Float64}, accuracy::Float64)
        self.L = M
        self.accuracy = accuracy
    end

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
                    yhat = vcat(yhat, sensor.measurement(x))
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
        function set0(M::Matrix{Float64}, accuracy::Float64)
            set(self, M, accuracy)
        end
        self.set = set0
        function system0()
            system(self)
        end
        self.system = system0
        return self
    end
end