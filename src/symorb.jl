module symorb
using GAP
using JSON
using OffsetArrays
using LinearAlgebra
using Plots
using Optim
using TimerOutputs

const to = TimerOutput()

include("globals.jl")
include("initdata.jl")
include("matrices.jl")
include("minimization.jl")
include("action.jl")
include("projectors.jl")


export exec


function exec(file)
    global Ω, K
    reset_timer!(to)

    read_init(file)

    A = (project ∘ random_starting_path)()

    results = path_minimize(A)
    show(to)
    cb(results)
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
    pl = plot(aspect_ratio=:equal)
    for i in 1:N
        plot!(pl, [[path[j][i][d] for j in 0:steps+1] for d in 1:dim]..., label="Body $i", linewidth=3,  aspect_ratio=:equal)
    end
    display(plot(pl))
    return false
end



end # module symorb
