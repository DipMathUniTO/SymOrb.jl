"""
    initialize(data::Dict)::Problem
    initialize(filename::String)::Problem

Initialize the problem with the data given in the dictionary `data` or in the file `filename`.
"""
function initialize(data::Dict)::Problem
    # Load GAP package
    Packages.load("$(@__DIR__)" * "/gap")
    
    # Set up the GAP variables
    GG.NOB = N = data["NOB"]
    GG.dim = dim =  data["dim"]
    at = data["action_type"]
    action_type = ActionType(at)

    GG.kern = evalstr(string(data["kern"]))

    GG.rotV = evalstr(data["rotV"])
    GG.rotS = evalstr(data["rotS"])
    if action_type != Cyclic
        GG.refV = evalstr(data["refV"])
        GG.refS = evalstr(data["refS"])
    else 
        GG.refV = evalstr( "[]")
        GG.refS = evalstr( "()")
    end
    # Create the GAP Symmetry Group
    GG.LSG = GG.LagSymmetryGroup(at, N, GG.kern, GG.rotV, GG.rotS, GG.refV, GG.refS)

    # Generate the elements of the symmetry group
    elements = g2j(GG.MinorbInitElements(GG.LSG), recursive=false)

    # Extract the elements of the symmetry group
    kerT = GroupElement.(elements[1])
    g = GroupElement.(elements[2])

    H0, H1 = if action_type != Cyclic
        GroupElement.((elements[3], elements[4]))
    else
        (GroupElement(), GroupElement())
    end

    m = convert.(Float64, data["m"])    # the masses
    F = data["F"]::Int                  # number of Fourier series terms
    dims = (F, N, dim)
    steps = if haskey(data, "steps")   # number of steps in the discretization of time [0,π]    
        data["steps"]::Int
    else
        2 * F
    end            

    Ω = Matrix{Float64}(hcat(data["Omega"]...)')         # generator of angular velocity
 
    # If the denominator of the potential energy is given, use it. Otherwise, use the identity function

    if haskey(data, "denominator")
        f = @RuntimeGeneratedFunction Meta.parse("x -> "*data["denominator"])
    else 
        f = identity
    end
    
    G = SymmetryGroup(action_type, kerT, g, H0, H1)
    R = nth_body(m, dims)
    Ri = remove_nth_body(dims)
    Π = Ri * project(G, dims) * R
    K = R' * kinetic_matrix(Ω, m, dims) * R      # Compute the kinetic energy matrix
    
    dx_dA = compute_dx_dA(F, steps)
    dA_dx = compute_dA_dx(F, steps)

    
    x_to_A = Ri * kron(dA_dx, I(N*dim)) 
    A_to_x = kron(dx_dA, I(N*dim)) * R
    
    # Matrix to integrate the potential
    Zu = compute_integration_factors(steps)
    Zg = Matrix((kron(Zu, ones(N*dim)) .* A_to_x)') 
    Zh = Matrix((kron(Zu', I(N*dim)) * A_to_x )')

    return Problem(N, dim, F, steps, G, m, f, K, A_to_x, x_to_A, Π, R, Ri, Zu, Zg, Zh, data)

end

initialize(filename::String)::Problem = initialize(TOML.parsefile(filename))

