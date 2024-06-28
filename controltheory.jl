module ControlTheory

    include("nodes.jl")
    include("observers.jl")
    include("controllers.jl")
    include("filters.jl")

    export Pump
    export MeasurementSensor
    export StateModel

    export LuenbergerObserver
    export optimize

    export ModelPredictiveController

    export ExponentialFilter

    using DifferentialEquations
    using Optimization
    using OptimizationNLopt

    function optimize(observer::OptimizableObserver, tSpan, u0, interval; maxtime=20)
        pF = [[],-1]
        function optFunction(du, u, p, t)
            yhat = nothing
            for i in 1:size(observer.inputs)[1]
                sensor = observer.inputs[i]  
                if (yhat == nothing)
                    yhat = sensor.measurement(u)
                else 
                    yhat = vcat(yhat, sensor.measurement(u))
                end
            end

            observer.model(du, u, p, t)

            mat = p[1] * (p[2] - yhat)
            i = 1
            while (i <= size(du)[1])
                du[i] = du[i] + mat[i]
                i = i + 1
            end
        end

        cbs = []
        tStops = []
        i = 1
        while i*interval <= tSpan[2]
            timeIndex= i*interval
            preTimeIndex = (i-1)*interval
            function c(u, t, intregrator)
                if (t == timeIndex)
                    println("Optimization at t=$(t)...")
                end
                return t==timeIndex
            end
            function a!(intregrator)
                ti = time()
                map = Dict()
                tv  = intregrator.u
                ym = nothing
                for i in 1:size(observer.inputs)[1]
                    sensor = observer.inputs[i]  
                    if (ym == nothing)
                        ym = sensor.measurement(tv)
                    else 
                        ym = vcat(ym, sensor.measurement(tv))
                    end
                end
                function opt(x, p)
                    tF = reshape(x, (size(tv)[1], size(ym)[1]))
                    optProblem = ODEProblem(optFunction, u0, (0,timeIndex))
                    eSol = solve(optProblem, Tsit5(), p=[tF, ym], verbose=false)
                    ev = eSol(timeIndex)
                    err = (tv-ev)/tv
                    sq = sum(e -> e^2, err)
                    if (pF[2] >= 0)
                        peSol = solve(optProblem, Tsit5(), p=[pF[1], ym], verbose=false)
                        pev = peSol(timeIndex)
                        peErr = (tv-pev)/tv
                        peSq =  sum(e -> e^2, peErr)
                        map[sq] = [sq, peSq]
                    else 
                        map[sq] = [sq, 0]
                    end
                    return sq
                end
                s = size(ym)[1] * size(tv)[1]
                bounds = [fill(0,s),fill(1,s)]
                x0 = bounds[1]
                p = [1.0]
                funca = OptimizationFunction(opt)
                probl = Optimization.OptimizationProblem(funca, x0, p, lb = bounds[1], ub = bounds[2])
                sol = solve(probl, NLopt.G_MLSL(), local_method=NLopt.LN_NELDERMEAD(), maxtime=maxtime)
                nF = reshape(sol.u, (size(intregrator.u)[1], size(ym)[1]))
                s = map[sol.objective]
                if (pF[2] < 0) 
                    pF[1] = reshape(sol.u, (size(intregrator.u)[1], size(ym)[1]))
                    pF[2] = s[1]
                else 
                    tw = s[1] + s[2]
                    itw = (tw/s[2]) + (tw/(s[1]))
                    v = ((tw/s[2])/itw)*pF[1] + ((tw/(s[1]))/itw)*(nF)
                    pF[1] = v
                    acc = (intregrator.t/tSpan[2])*pF[2] + (interval/tSpan[2])*(s[1])
                    pF[2] = acc
                end
                println("Optimization: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
                observer.set(pF[1], pF[2])

                if ((timeIndex + interval) > tSpan[2])
                    println("Optimization: Finished with $(pF[1]) with accuracy $(pF[2]).")
                    observer.set(pF[1], pF[2])
                end
            end
            cb = DiscreteCallback(c,a!) 
            push!(cbs,cb)
            push!(tStops, timeIndex)
            i = i + 1
        end

        function trueModel(du, u, p, t)
            for i in 1:size(observer.inputs)[1]
                sensor = observer.inputs[i]  
                sensor.parent.model(du,u,p,t)
            end
        end


        return (trueModel, CallbackSet(cbs...), tStops)
    end
end;