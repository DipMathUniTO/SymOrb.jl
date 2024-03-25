function path_minimize(Γ)

    v = flatten(Γ)
    show_infos = false
    
    @timeit to "Optimize" begin
        res = optimize(v_action,  v_∇action,  v, BFGS(), Optim.Options(iterations=1000, show_trace=false, extended_trace=true,
            callback=x -> begin
            if (show_infos)
            if (x.iteration % 10 == 0)
                println(x.iteration, "\t", x.value, "\t", x.g_norm)
                cb(x.metadata["x"])
                end
            end
                return false
            end); inplace = false)
    end
    return res.minimizer

end