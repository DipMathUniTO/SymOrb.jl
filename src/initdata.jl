
GG = GAP.Globals
g2j = GAP.gap_to_julia


function E(n)
   return real(exp(2*π * im / n))
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

    # load GAP package
    GAP.Packages.load("$(@__DIR__)" * "/gap")


    global N = data["NOB"]
    GG.NOB = N
    
    global dim = data["dim"]
    GG.dim = dim
    
    at = data["action_type"]
    global action_type = ActionType(at)

    GG.kern = GAP.evalstr(data["kern"])
    GG.rotV = GAP.evalstr(data["rotV"])
    GG.rotS = GAP.evalstr(data["rotS"])
    GG.refV = GAP.evalstr(data["refV"])
    GG.refS = GAP.evalstr(data["refS"])

    GG.LSG = GG.LagSymmetryGroup(at, N, GG.kern, GG.rotV, GG.rotS, GG.refV, GG.refS)

    minorb_elements = g2j(GG.MinorbInitElements(GG.LSG), recursive=false)

    global is_typeR =  if dim==3 GG.IsTypeR(GG.LSG) else false end
    global kerT = perm_from_gap(minorb_elements[1])
    global g = perm_from_gap(minorb_elements[2])

    if action_type != Cyclic
        global H_0 = perm_from_gap([minorb_elements[3]])[1]
        global H_1 = perm_from_gap([minorb_elements[4]])[1]
        global H_0 = perm_from_gap([minorb_elements[3]])[1]
        global H_1 = perm_from_gap([minorb_elements[4]])[1]
    else
        global H_0 = G()
        global H_1 = G()
        global H_0 = G()
        global H_1 = G()
    end

    global cyclic_order = length(g)
    global cyclic_order = length(g)

    global m = data["m"]           # the masses
    global F = data["F"]           # number of Fourier series terms
    global steps = 2 * F           # number of steps in the discretization of time [0,1]
    global dt = π / (steps+1)      # time step
    global Ω = hcat(data["Omega"]...)

    global Ω2 = Ω * Ω
    global dx_dAk = compute_dx_dAk()
    global Ω2 = Ω * Ω
    global dx_dAk = compute_dx_dAk()
    
    global Id = [if i == j m[i] * I(dim) else zeros(dim, dim) end for i in 1:N, j in 1:N]
    global Id = [if i == j m[i] * I(dim) else zeros(dim, dim) end for i in 1:N, j in 1:N]

    global K = K_linear()
    global K = K_linear()
    
    if (!iszero(Ω))
        K +=  Ω2 * K_centrifugal() + Ω * K_coriolis()
    end

    return 
end

