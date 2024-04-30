
GG = GAP.Globals
g2j = GAP.gap_to_julia


function perm_from_gap(mat)
    gs = Vector{G}()
    for el in mat
        GG.tuple = GapObj(el)
        perm::Vector =  g2j(@gap Permuted([1 .. NOB], tuple[2]^(-1)))
        matrix::Matrix = hcat(g2j(@gap tuple[1])...)
        push!(gs, G(perm, matrix))
    end
    return gs
end


function LGS_from_config(data::AbstractDict)
    global kerT, g, H_0, H_1, m, F, Ω, Ω2, α, steps, N, dim, action_type, cyclic_order, is_typeR, isΩ, K

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

    rotV = GapObj(data["rotV"], recursive=true)
    rotS = GAP.evalstr(data["rotS"])
    refV = GapObj(data["refV"], recursive=true)
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
    α = 1.                  # Kepler exponent
    steps = 2 * F           # number of steps in the discretization of time [0,1]
    Ω = hcat(data["Omega"]...)
    isΩ = !iszero(Ω)
    Ω2 = Ω * Ω
    K = K_linear()
    if (isΩ)
        K +=  spatial_mult(Ω2, K_centrifugal()) + spatial_mult(Ω, K_coriolis())
    end

end

