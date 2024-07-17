struct Configuration{domain, N, V<:AbstractArray}
    parent::V
    offsets::NTuple{N, Int}
end

function Configuration{domain, N}(v::V, first_indices=ntuple(Returns(1), N)) where {domain, N, V<:AbstractArray{<:Any, N}}
    offsets = first_indices .- ntuple(Returns(1), N)
    return Configuration{domain, N, V}(v, offsets)
end

const Path = Configuration{:time, 3}
const Coefficients = Configuration{:fourier,3}
const Config = Configuration{:space, 2}

function Base.getindex(c::Configuration, is...)
    I = is .- c.offsets[1:length(is)]
    # TODO: handle offset arrays here
    colons = ntuple(Returns(:), ndims(c.parent) - length(I))
    return view(c.parent, colons..., I...)
end

function Base.setindex!(c::Configuration, val, is...)
    c[is...] .= val
    return c
end

Base.parent(c::Configuration) = c.parent
Base.size(c::Configuration) = size(parent(c))


# TODO: add helpers as needed

# sum(c[1:5], dims=2)
# @. c[3] *= 2 (or overload `Base.setindex!` but test `@allocated`)
# TensorCast or TensorOperations (inplace operations on views)
# for el in c[] (iteration)

# TODO port emboss & co. to new format

# to_vector(c::Configuration) = vec(parent(c))

# function to_matrix(c::Configuration)
#     m = permutedims(parent(c), (1, 3, 5, 2, 4, 6))
#     sz = size(m)
#     return reshape(m, prod(sz[1:3]), prod(sz[4:6]))
# end
