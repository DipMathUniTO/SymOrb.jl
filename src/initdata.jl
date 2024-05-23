
GG = GAP.Globals
g2j = GAP.gap_to_julia


function E(n)
   return real(exp(2π * im * n / N))
end

function perm_from_gap(list)
    gs = Vector{G}()
    for el ∈ list
        GG.tuple = GapObj(el)
        perm::Vector =  g2j(@gap Permuted([1 .. NOB], tuple[2]^(-1)))
        matrix::Matrix = eval.(Meta.parse.(hcat(g2j(@gap ConvertToString(tuple[1]))...)))
        
        push!(gs, G(perm, matrix))
    end
    return gs
end


function LSG_from_config(data::AbstractDict)
    global kerT, g, H_0, H_1, m, F, Ω, Ω2, dt, Id, steps, N, dim, action_type, cyclic_order, is_typeR, K, dx_dAk

    # load GAP package
    GAP.Packages.load("$(@__DIR__)" * "/gap")


    N = data["NOB"]
    dim = data["dim"]
    at = data["action_type"]

    action_type = ActionType(at)
    if data["kern"] == "TrivialKerTau"
        kern = GG.TrivialKerTau(dim)
    else
        kern = GapObj(data["kern"], recursive=true)
    end

    rotV = GAP.evalstr(data["rotV"])
    rotS = GAP.evalstr(data["rotS"])
    refV = GAP.evalstr(data["refV"])
    refS = GAP.evalstr(data["refS"])

    
    GG.dim, GG.NOB, GG.kern, GG.rotV, GG.rotS, GG.refV, GG.refS = dim, N, kern, rotV, rotS, refV, refS
    
    GG.LSG = GG.LagSymmetryGroup(at, N, kern, rotV, rotS, refV, refS)
    minorb_elements = g2j(GG.MinorbInitElements(GG.LSG), recursive=false)
    is_typeR =  if dim==3 GG.IsTypeR(GG.LSG) else false end

    kerT = perm_from_gap(minorb_elements[1])
    g = perm_from_gap(minorb_elements[2])

    if action_type != Cyclic
        H_0 = perm_from_gap([minorb_elements[3]])[1]
        H_1 = perm_from_gap([minorb_elements[4]])[1]
    else
        H_0, H_1 = G(), G()
    end

    cyclic_order = length(g)

    m = data["m"]              # the masses
    F = data["F"]              # number of Fourier series terms
    steps = 2 * F           # number of steps in the discretization of time [0,1]
    dt = π / (steps+1)      # time step
    Ω = hcat(data["Omega"]...)

    Ω2 = Ω * Ω
    dx_dAk = compute_dx_dAk()
    
    Id = [if i == j I(dim) else zeros(dim, dim) end for i in 1:N, j in 1:N]

    K = K_linear()
    
    if (!iszero(Ω))
        K +=  Ω2 * K_centrifugal() + Ω * K_coriolis()
    end

    return 
end

