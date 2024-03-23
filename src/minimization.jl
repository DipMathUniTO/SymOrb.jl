function path_minimize(Γ)

    v = flatten(Γ)

    @timeit to "Optimize" begin
        res = optimize(v_action, v, GradientDescent(), Optim.Options(iterations=1000, show_trace=false, extended_trace=true,
            callback=x -> begin
                println(x.iteration, "\t", x.value, "\t", x.g_norm)

                if (x.iteration % 2 == 0)
                    cb(x.metadata["x"])
                end

                return false
            end))
    end
    return res.minimizer

end