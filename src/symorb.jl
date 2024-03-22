module symorb
using GAP
using JSON
using OffsetArrays
using LinearAlgebra
using Plots
using ForwardDiff
using Calculus
using Optim
using ForwardDiff



include("globals.jl")
include("initdata.jl")

include("action.jl")
include("projectors.jl")



export exec


function exec(file)
    global Ω

    read_init(file)
    A = random_starting_path()
    A = project(A)
    results = path_minimize(A)
    @show results
    return
end

function random_starting_path()
    A = [[rand(Float64, dim) .- 0.5 for i in 1:N] for j in 0:F+1]
    return OffsetArrays.Origin(0)(A)
end


function circular_starting_path()
    A::Path = OffsetArrays.Origin(0)([[zeros(dim) for i in 1:N] for j in 1:steps+2])

    for h in 0:steps+1
        for i in 1:N
            A[h][i][1] = cos(π*h/(steps+1) + (i-1)*2*π/N)
            A[h][i][2] = sin(π*h/(steps+1) + (i-1)*2*π/N)
        end
    end
    return OffsetArrays.Origin(0)(inverse_fourier(A))
end


function cb(v)
    path = build_path((project ∘ emboss)(v))
    pl = plot(
        [path[i][1][1] for i in 0:steps + 1],
        [path[i][1][2] for i in 0:steps+1], label="Body 1", markershape=:circle)
    plot!(pl, [path[i][2][1] for i in 0:steps+1],
        [path[i][2][2] for i in 0:steps + 1], label="Body 2", markershape= :circle)
    plot!(pl, [path[i][3][1] for i in 0:steps + 1],
        [path[i][3][2] for i in 0:steps + 1], label="Body 3", markershape= :circle)
    display(plot(pl))
    return false
end


function custom_gradient_descent(Γ)
     v = flatten(Γ)
     p = 1
     grad = Calculus.gradient(v_action, v)
     grad = (flatten∘project∘emboss)(grad)
     @show grad

    nm = norm(grad)

    println(nm)
    old = v_action(v)
    step = p
    iter = 0

    cb(v)
    while true
        v_new = (flatten ∘ project ∘ emboss)(v - step * grad)

        new = v_action(v_new)
        if(new > old)
            step /= 2
            continue
        end

        grad_new = (flatten ∘ project ∘ emboss)(Calculus.gradient(v_action, v_new))

        if(iter % 5 == 0) cb(v_new) end
        iter +=1
        old = new
        v = v_new
        grad = grad_new
        nm = norm(grad)
        println(step, " ",new,"      ", nm)

        step = min(p, step * 2)
        iter += 1
        if nm < 1e-6
            break
        end
        if iter > 1e6
            break
        end
        if step < 1e-6
            break
        end
    end
end


function path_minimize(Γ)

    v = flatten(Γ)


        res = optimize(v_action, v, GradientDescent(), Optim.Options(iterations=1000, show_trace=false, extended_trace=true, callback=x -> begin
            println(x.iteration, "\t", x.value, "\t", x.g_norm)

            if (x.iteration % 10 == 0) cb(x.metadata["x"]) end

            return false
        end))
        v = res.minimizer #(flatten ∘ project ∘ emboss)(res.minimizer)

        #cb(v)

    #x, info = prima(v_action, v; iprint=PRIMA.MSG_FEVL, nonlinear_eq = constraint)
end


end # module symorb
