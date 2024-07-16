# TensorOperations
struct Configuration{domain, N, V<:AbstractArray}
    parent::V
    # TODO: add additional fields
end

function Configuration{domain, N}(v::V) where {domain, N, V<:AbstractArray{<:Any, N}}
    return Configuration{domain, N, V}(v)
end

const Path = Configuration{:time, 3}
const Coefficients = Configuration{:frequency, 2}

# f(p::Path)

# path::Configuration
# freqs::Configuration

function Base.getindex(c::Configuration, is...)
    I = reverse(is)
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
