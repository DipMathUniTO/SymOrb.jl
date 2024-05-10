#  action of the inverse of g ∈ G on v ∈ ℝᵈᴺ
ϕg_inv(v::Config)::Config =
    if (cyclic_order > 1)
        ϕ(g[1], v)
    else
        v
    end

# action of g ∈ G on v ∈ ℝᵈᴺ
ϕg(v::Config)::Config =
    if (cyclic_order > 1)
        ϕ(g[cyclic_order-1], v)
    else
        v
    end

# Projection onto the center of mass
π_c(v::Config)::Config =  begin
    xc = sum(m[i] * v[i] for i ∈ axes(m, 1)) / sum(m)

    return [v[i] - xc for i ∈ axes(m,1)]
end


# action of h ∈ H on v ∈ ℝᵈ
ϕ(h::G, v::Config)::Config =
    [h.M * v[h.σ[i]] for i ∈ axes(h.σ, 1)]

# Projection onto the subspace spanned by the action of H on v
π_H(H::Vector{G}, v::Config)::Config =
    sum(ϕ(H[i], v) for i ∈ axes(H, 1)) / length(H)

# Projection onto ker_T
π_kerT(v::Config)::Config = π_H(kerT, v)

# Simplified formulas since H_0 and H_1 are cyclic of order 1
π_01(v::Config, w::Config)::Tuple{Config,Config} = ((v + ϕ(H_0, v))/2,    (w + ϕ(H_1, w)) / 2)

# Case of cyclic action
π_g(v::Config, w::Config)::Tuple{Config, Config} =  ((v + ϕg_inv(w)) / 2,  (w + ϕg(v)) / 2)


π_a(v::Config, w::Config)::Tuple{Config,Config} =
    if (action_type == Cyclic)
        π_g( π_c(v), π_c(w) )
    else
        π_01( π_c(v), π_c(w) )
    end

function project(A::Coefficients)::Coefficients
    (A[0], A[F+1]) = π_a(A[0], A[F+1])
    A = π_kerT.(π_c.(A))
    return A
end

