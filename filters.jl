abstract type StateFilter <: StateProcess end


mutable struct ExponentialFilter <: StateFilter


    A::Matrix{Float64}

    parent::StateProcess
    T::Float64
    τ::Float64

    system::Function
    model::Function

    function system(self::ExponentialFilter)
        function system0(du, u, p, t)
            dy = du * self.A
            y = u * self.A
            B = transpose(self.A)

            self.parent.system()(du,u,p,t)
            dx = du * self.parent.A
            x = u * self.parent.A
            α = exp((-self.T)/self.τ)
            dy .=  α*dy + (1-α)*dx

            u .= y*B - (u*self.A)*B + u
            du .= dy*B - (du*self.A)*B + du
        end
    end

    function ExponentialFilter(A::Matrix{Float64}, parent::StateProcess, T::Float64, τ::Float64)
        self = new(A, parent, T, τ)
        function system0()
            system(self)
        end
        self.system = system0
        self.model = parent.model
        return self
    end
end