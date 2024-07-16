function E(n)
   return real(exp(2*π * im / n))
end

"""
    to_julia(el)

Create a GroupElement from the GAP result.
"""
function to_julia(el)::GroupElement
    GG.tuple = GapObj(el)
    perm::Permutation =  g2j(@gap Permuted([1 .. NOB], tuple[2]^(-1)))
    matrix::Rotation = eval.(Meta.parse.(hcat(g2j(@gap ConvertToString(tuple[1]))...)))  
    return GroupElement(perm, matrix)    
end

"""
    initialize(data::T) where {T<:AbstractDict}

Initialize the problem with the data given in the dictionary.
"""
function initialize(data::T) where {T<:AbstractDict}
    # Load GAP package
    Packages.load("$(@__DIR__)" * "/gap")
    
    # Set up the GAP variables
    GG.NOB = N = data["NOB"]
    GG.dim = dim =  data["dim"]
    at = data["action_type"]
    action_type = ActionType(at)
    GG.kern = evalstr(data["kern"])
    GG.rotV = evalstr(data["rotV"])
    GG.rotS = evalstr(data["rotS"])
    GG.refV = evalstr(data["refV"])
    GG.refS = evalstr(data["refS"])

    # Create the GAP Symmetry Group
    GG.LSG = GG.LagSymmetryGroup(at, N, GG.kern, GG.rotV, GG.rotS, GG.refV, GG.refS)

    # Generate the elements of the symmetry group
    elements = g2j(GG.MinorbInitElements(GG.LSG), recursive=false)

    # Extract the elements of the symmetry group
    kerT = to_julia.(elements[1])
    g = to_julia.(elements[2])
    H0, H1 = if action_type != Cyclic
        to_julia.((elements[3], elements[4]))
    else
        (GroupElement(), GroupElement())
    end

    m = convert.(Float64, data["m"])    # the masses
    F = data["F"]::Int64                # number of Fourier series terms
    steps = 2 * F                       # number of steps in the discretization of time [0,1]    
    dx_dAk = compute_dx_dAk(F, steps)   # the derivative of the path with respect to the Fourier coefficients
    M = [if i == j m[i] * I(dim) else zeros(dim, dim) end for i in 1:N, j in 1:N] # Masses matrix
    
    # Compute the kinetic energy matrix
    Ω = hcat(data["Omega"]...)'         # generator of angular velocity
    K = K_linear(N, F, dim, M)          # linear part of the kinetic energy matrix 
    if (!iszero(Ω))
        K +=  K_centrifugal(Ω*Ω, N, F, dim, M) + K_coriolis(Ω, N, F, dim, M)
    end

    # If the denominator of the potential energy is given, use it. Otherwise, use the identity function
    f = if haskey(data, "denominator")
            eval(Meta.parse("x -> " * data["denominator"]))
        else 
            x -> x
        end

    G = SymmetryGroup(action_type, kerT, g, H0, H1)
    return Problem( (F=F, N=N, dim=dim), G, m, f, K, dx_dAk)

end


initialize(file::String) = initialize(parsefile(file))
