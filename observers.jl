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
        function system0()
            system(self)
        end
        self.system = system0
        return self
    end
end



# Optimization
# abstract type SE_Optimization end

# mutable struct SE_NelderMead <: SE_Optimization

#     interval::Int
#     maxTime::Float64

#     optimize::Function

#     function optimize(self::SE_NelderMead, observer::OptimizableObserver, u0, tSpan)
#         function format(observer::LuenbergerObserver, stateEstimator::StateEstimator, x::Array)
#             testArr = stateEstimator.measurement(stateEstimator.trueProblem.u0)
#             formatted = reshape(x, (size(stateEstimator.trueProblem.u0)[1], size(testArr)[1]))
#             return formatted
#         end
        
#         pF = [[],-1]
#         function optFunction(du, u, p, t)
#             stateEstimator.estimatorProblem.f.f(du,u,p,t)

#             x = u * observer.A
#             dx = du * observer.A
#             B = transpose(self.A)
#             yhat = nothing
#             for i in 1:size(observer.inputs)[1]
#                 sensor = observer.inputs[i]  
#                 if (yhat == nothing)
#                     yhat = sensor.measurement(x)
#                 else 
#                     yhat = vcat(yhat, sensor.measurement(x))
#                 end
#             end

#             observer.model(dx, x, p, t)

#             mat = p[1] * (ym - yhat)
#             i = 1
#             while (i <= size(dx)[1])
#                 dx[i] = dx[i] + mat[i]
#                 i = i + 1
#             end

#             u .= x*B - (u*self.A)*B + u
#             du .= dx*B - (du*self.A)*B + du
#         end
#         optProblem = ODEProblem(optFunction, u0, tSpan)

#         cbs = []
#         tStops = []
#         i = 1
#         while i*self.interval <= tSpan[2]
#             timeIndex = i*self.interval
#             preTimeIndex = (i-1)*self.interval
#             function c(u,t,intregrator) 
#                 if (t == timeIndex)
#                     println("Optimizing at $(timeIndex)...")
#                 end
#                 t == timeIndex
#             end
#             function a!(intregrator)
#                 ti = time()
#                 map = Dict()
#                 function opt(x, p) 
#                     tv  = intregrator.u
#                     ym = nothing
#                     for i in 1:size(observer.inputs)[1]
#                         sensor = observer.inputs[i]  
#                         if (ym == nothing)
#                             ym = sensor.measurement(intregrator.u)
#                         else 
#                             ym = vcat(ym, sensor.measurement(intregrator.u))
#                         end
#                     end
#                     tF = reshape(x, (size(intregrator.u)[1], size(ym)[1]))
#                     eSol = solve(optProblem, Tsit5(), p=[tF, ym], verbose=false)
#                     ev = eSol(timeIndex)
#                     err = tv-ev
#                     sq = sum(e -> e^2, err)
#                     if (pF[2] >= 0)
#                         peSol = solve(optProblem, Tsit5(), p=[pF[1], ym], verbose=false)
#                         pev = peSol(timeIndex)
#                         peErr = tv-pev
#                         peSq =  sum(e -> e^2, peErr)
#                         map[sq] = [sq, peSq]
#                     else 
#                         map[sq] = [sq, 0]
#                     end
#                     return sq
#                 end
#                 s = size(ym)[1] * size(intregrator.u)[1]
#                 return [fill(0,s),fill(1,s)]
#                 x0 = bounds[1]
#                 p = [1.0]
#                 funca = OptimizationFunction(opt)
#                 probl = Optimization.OptimizationProblem(funca, x0, p, lb = bounds[1], ub = bounds[2])
#                 sol = solve(probl, NLopt.LN_NELDERMEAD(), maxtime=self.maxTime)
#                 nF = reshape(sol.u, (size(intregrator.u)[1], size(ym)[1]))
#                 s = map[sol.objective]
#                 if (pF[2] < 0) 
#                     pF[1] = reshape(sol.u, (size(intregrator.u)[1], size(ym)[1]))
#                     pF[2] = s[1]
#                 else 
#                     tw = s[1] + s[2]
#                     itw = (tw/s[2]) + (tw/(s[1]))
#                     v = ((tw/s[2])/itw)*pF[1] + ((tw/(s[1]))/itw)*(nF)
#                     pF[1] = v
#                     acc = ((tw/s[2])/itw)*pF[2] + ((tw/(s[1]))/itw)*(s[1])
#                     pF[2] = acc
#                 end
#                 println("SEOptimization: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
#             end
#             cb = DiscreteCallback(c,a!) 
#             push!(cbs,cb)
#             push!(tStops, timeIndex)
#             i = i + 1
#         end

#         function trueModel(du, u, p, t)
#             tP = []
#             for i in 1:size(observer.inputs)[1]
#                 sensor = observer.inputs[i]  
#                 if (!(sensor.parent in tP))
#                     push!(tP, sensor.parent)
#                 end
#             end
#             for i in 1:size(tP)
#                 tP[i].system()(du,u,p,t)
#             end
#         end

#         trueSolution = solve(trueModel, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

#         println("SEOptimization: Finished with $(pF[1]) with accuracy $(pF[2]).")
#         observer.L = pF[1]
#         observer.accuracy = pF[2]
#         return pF
#     end

#     function SE_NelderMead(interval; maxTime=20)
#         self = new(interval, maxTime)
#         function optimize0(observer::OptimizableObserver)
#             optimize(self, observe)
#         end
#         self.optimize = optimize0
#         return self
#     end
# end

# mutable struct SE_PolyChaos <: SE_Optimization

#     interval::Int
#     n::Int
#     d::Int

#     optimize::Function

#     function optimize(self::SE_PolyChaos, stateEstimator::StateEstimator)
#         pF = [[],-1]
#         function optFunction(du, u, p, t)
#             stateEstimator.estimatorProblem.f.f(du,u,p,t)

#             y_hat = stateEstimator.measurement(u)
#             stateEstimator.observer.set(stateEstimator, p[1], 0.0)
#             stateEstimator.observer.observe(stateEstimator, p[2], y_hat, du, u, p, t)
#         end
#         optProblem = ODEProblem(optFunction, stateEstimator.estimatorProblem.u0, stateEstimator.tSpan)

#         cbs = []
#         tStops = []
#         i = 1
#         while i*self.interval <= stateEstimator.tSpan[2]
#             timeIndex = i*self.interval
#             preTimeIndex = (i-1)*self.interval
#             function c(u,t,intregrator) 
#                 if (t == timeIndex)
#                     println("Optimizing at $(timeIndex)...")
#                 end
#                 t == timeIndex
#             end
#             function a!(intregrator)
#                 ti = time()
#                 map = Dict()
#                 function opt(x) 
#                     ym = stateEstimator.measurement(intregrator.u)
#                     tF = stateEstimator.observer.format(stateEstimator, [x...])
#                     tv  = intregrator.u
#                     eSol = solve(optProblem, Tsit5(), p=[tF, ym], maxiters=1e4, verbose=false)
#                     ev = eSol(timeIndex)
#                     err = tv-ev
#                     sq = sum(e -> e^2, err)
#                     if (pF[2] >= 0)
#                         peSol = solve(optProblem, Tsit5(), p=[pF[1], ym], maxiters=1e4, verbose=false)
#                         pev = peSol(timeIndex)
#                         peErr = tv-pev
#                         peSq =  sum(e -> e^2, peErr)
#                         map[sq] = [sq, peSq]
#                     else 
#                         map[sq] = [sq, 0]
#                     end
#                     return sq
#                 end
#                 bounds = stateEstimator.observer.bounds(stateEstimator)
#                 uniform01 = Uniform01OrthoPoly(self.d)
#                 samples = []
#                 j = 1
#                 while j <= size(bounds[1])[1]
#                     push!(samples, evaluatePCE(convert2affinePCE(bounds[1][j],bounds[2][j],uniform01), sampleMeasure(self.n,uniform01),uniform01))
#                     j = j + 1
#                 end
#                 lErr = Inf
#                 lEle = []
#                 z = 1
#                 for elements in Iterators.product(samples...)
#                     sq = opt(elements)
#                     if (sq < lErr) 
#                         lErr = sq
#                         lEle = elements
#                     end
#                     z = z + 1
#                 end
#                 nF = stateEstimator.observer.format(stateEstimator, [lEle...])
#                 s = map[lErr]
#                 if (pF[2] < 0) 
#                     pF[1] = stateEstimator.observer.format(stateEstimator, [lEle...])
#                     pF[2] = s[1]
#                 else 
#                     tw = s[1] + s[2]
#                     itw = (tw/s[2]) + (tw/(s[1]))
#                     v = ((tw/s[2])/itw)*pF[1] + ((tw/(s[1]))/itw)*(nF)
#                     pF[1] = v
#                     acc = ((tw/s[2])/itw)*pF[2] + ((tw/(s[1]))/itw)*(s[1])
#                     pF[2] = acc
#                 end
#                 println("SEOptimization: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
#             end
#             cb = DiscreteCallback(c,a!) 
#             push!(cbs,cb)
#             push!(tStops, timeIndex)
#             i = i + 1
#         end

#         trueSolution = solve(stateEstimator.trueProblem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

#         println("SEOptimization: Finished with $(pF[1]) with accuracy $(pF[2]).")
#         observer.L = pF[1]
#         observer.accuracy = pF[2]
#         return pF
#     end

#     function SE_PolyChaos(interval, n; d=6)
#         self = new(interval, n, d)
#         function optimize0(stateEstimator::StateEstimator)
#             optimize(self, stateEstimator)
#         end
#         self.optimize = optimize0
#         return self
#     end
# end