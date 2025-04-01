module BracketAlgebras

using AbstractAlgebra
# import Base: show, +, -, *, ^, ==, <, >, inv, isone, iszero, one, zero, rand, deepcopy_internal, hash
using Bijections
using Combinatorics

abstract type AbstractBracketAlgebra{T} <: Ring end
abstract type AbstractBracketAlgebraElem{T} <: RingElem end

import AbstractAlgebra.gens

export AbstractBracketAlgebra, BracketAlgebra, d, n, Tabloid, is_standard, bracket_monomial, AbstractBracketAlgebraElem, BracketAlgebraElem, point_ordering, point_ordering!, point_labels, point_labels!, standard_violation, straightening_sizyge, straighten, atomic_extensors, is_multilinear


# aux function for the lexicographic order of vectors extended from point_ordering
# point_ordering[1] < point_ordering[2] < ...
# See Sturmfels 2008 p.81
function _lt(point_ordering::Vector{Int})
    function lt(a::Integer, b::Integer)
        indexin(a, point_ordering)[1] < indexin(b, point_ordering)[1]
    end
    function lt(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
        return indexin(sort(a, lt=lt), point_ordering) < indexin(sort(b, lt=lt), point_ordering)
    end
    return lt
end

"""
     BracketAlgebra{T<:RingElement} <: AbstractBracketAlgebra{T}
 
 Bracket algebra with coefficients of type `T`.
 
 # Fields
 - `d::Int`: dimension of the underlying projective space P^`d`
 - `n::Int`: number of points considered
 - `R::MPolyRing{T}`: AbstractAlgebra multivariate polynomial ring representing the elements of the bracket algebra (for internal use)
 - `point_ordering::Vector{Int}`: ordering of the `n` point indices, that induces the tableaux order. `point_ordering[1] < point_ordering[2] < ... < point_ordering[n]`. 
 - `variables::Bijection{Vector{Int},<:MPolyRingElem{T}}`: bijection from brackets to the corresponding monomials in the polynomial ring `R`
 - `point_labels::Vector`: point labels for each point
 """
mutable struct BracketAlgebra{T<:RingElement} <: AbstractBracketAlgebra{T}
    d::Int # dimension of the underlying projective space $\mathbb{P}^d$
    n::Int # number of points considered
    R::MPolyRing{T}
    point_ordering::Vector{Int} # ordering of the n pointsjulia, that induces the tableaux order. point_ordering[1] < point_ordering[2] < ... < point_ordering[n]. See Sturmfels 2008 p.81
    variables::Bijection{Vector{Int},<:MPolyRingElem{T}}
    point_labels::Vector  # point labels for each point
end

"""
    BracketAlgebra(n::Integer, d::Integer, point_ordering::AbstractVector{<:Integer}=collect(1:n), coefficient_ring::AbstractAlgebra.Ring=AbstractAlgebra.ZZ; point_labels::AbstractVector=collect(1:n))


    BracketAlgebra(n::Integer, d::Integer, point_ordering::AbstractVector{<:Integer}=collect(1:n), T::Type=Int; point_labels::AbstractVector=collect(1:n))

Construct a bracket algebra over the projective space P^d on n points with the given point ordering and point labels. If point_labels is a vector of integers, it is expected to be equal to `collect(1:n)`. `coefficient_ring` is the ring, which contains the coefficients of the bracket algebra.

# Examples
```jldoctest
julia> B = BracketAlgebra(4, 2, [1,2,3,4], point_labels = ["a","b","c","d"])
Bracket algebra over P^2 on 4 points with point ordering a < b < c < d and coefficient ring Integers.
```

```jldoctest
julia> B = BracketAlgebra(4, 2, [2,1,3,4], point_labels = ["a","b","c","d"])
Bracket algebra over P^2 on 4 points with point ordering b < a < c < d and coefficient ring Integers.
```

```jldoctest
julia> B = BracketAlgebra(4, 2, [2,1,3,4], AbstractAlgebra.GF(13), point_labels = ["a","b","c","d"])
Bracket algebra over P^2 on 4 points with point ordering b < a < c < d and coefficient ring Finite field F_13.
```
"""
function BracketAlgebra(n::Integer, d::Integer, point_ordering::AbstractVector{<:Integer}=collect(1:n), coefficient_ring::AbstractAlgebra.Ring=AbstractAlgebra.ZZ; point_labels::AbstractVector=collect(1:n))
    if sort(point_ordering) != collect(1:n)
        error("ordering needs to order the points 1:$n, but got $point_ordering")
    end
    if length(point_labels) != n
        error("point_labels must be of length $n")
    end
    # Add check for point_labels when it's a Vector of Integers
    if eltype(point_labels) <: Integer && point_labels != collect(1:n)
        error("When point_labels is a Vector of Integers, it must equal collect(1:n)")
    end

    brackets = sort(sort.(collect(combinations(1:n, d + 1)), lt=_lt(point_ordering)), lt=_lt(point_ordering)) # brackets are in the form [a,b,c,...] with a < b < c according to point_ordering. They are sorted in decreasing order wrt the tableaux order extrapolated from point_ordering.
    vars = AbstractAlgebra.variable_names("x#" => brackets)

    # the monomial ordering is extrapolated from vars[1] > vars[2] > vars[3]... As the brackets are sorted in increasing order wrt the tableaux order, the variables need to be reversed.
    R, x = polynomial_ring(coefficient_ring, reverse(vars), internal_ordering=:degrevlex)
    reverse!(x) # now x is in the same order as brackets

    variables = Bijection{Vector{Int},typeof(x[1])}()

    for (i, bracket) in enumerate(brackets)
        variables[bracket] = x[i]
    end

    # ordering = Groebner.DegRevLex(x)

    return BracketAlgebra{elem_type(coefficient_ring)}(d, n, R, point_ordering, variables, point_labels)
end

# extension of AbstractAlgebra.gens
# """
#     gens(B::BracketAlgebra)

# Return the generators of the bracket algebra `B` as elements of `B`. The list is sorted decreasing wrt the tableaux order and is thus consistent with the ordering of the variables in the polynomial ring `B.R`.

# # Examples
# ```jldoctest
# julia> gens(BracketAlgebra(4,2,[1,2,3,4]))
# 4-element Vector{BracketAlgebraElem{BigInt}}:
#  [2, 3, 4]
#  [1, 3, 4]
#  [1, 2, 4]
#  [1, 2, 3]
# ```

# ```jldoctest
# julia> gens(BracketAlgebra(4,2,[2,1,3,4]))
# 4-element Vector{BracketAlgebraElem{BigInt}}:
#  [1, 3, 4]
#  [2, 3, 4]
#  [2, 1, 4]
#  [2, 1, 3]
# ```
# """
gens(B::BracketAlgebra) = [B(B.variables(x)) for x in gens(B.R)]


"""
    point_labels(B::BracketAlgebra)

Return the point labels of the bracket algebra `B`.
"""
point_labels(B::BracketAlgebra) = B.point_labels

"""
    point_labels!(B::BracketAlgebra, new_labels::AbstractVector)

Update the point labels of the bracket algebra `B` and return the resulting bracket algebra. If `new_labels` is a vector of integers, it is expected to be equal to `collect(1:n(B))`.

# Examples
```jldoctest
julia> point_labels!(BracketAlgebra(6,3), ["a", "b", "c", "d", "e", "f"])
Bracket algebra over P^3 on 6 points with point ordering a < b < c < d < e < f and coefficient ring Integers.
```

```jldoctest
julia> point_labels!(BracketAlgebra(6,3), [2,1,3,4,5,"a"])
Bracket algebra over P^3 on 6 points with point ordering 2 < 1 < 3 < 4 < 5 < a and coefficient ring Integers.
```

```jldoctest
julia> point_labels!(BracketAlgebra(6,3), [2,1,3,4,5,6])
ERROR: When point_labels is a Vector of Integers, it must equal collect(1:n)
```
"""
function point_labels!(B::BracketAlgebra, new_labels::AbstractVector)
    if length(new_labels) != B.n
        error("new_labels must be of length $(B.n)")
    end

    if eltype(new_labels) <: Integer && new_labels != collect(1:n(B))
        error("When point_labels is a Vector of Integers, it must equal collect(1:n)")
    end

    B.point_labels = new_labels
    return B
end

"""
    d(B::BracketAlgebra)

Return the dimension of the projective space P^`d` underlying the bracket algebra `B`.

# Examples
```jldoctest
julia> d(BracketAlgebra(6,3))
3
```
"""
d(B::BracketAlgebra) = B.d

"""
    n(B::BracketAlgebra)

Return the number of points considered in the bracket algebra `B`.
```jldoctest
julia> n(BracketAlgebra(6,3))
6
```
"""
n(B::BracketAlgebra) = B.n

"""
    point_ordering(B::BracketAlgebra)

Return the ordering of the n points in the bracket algebra `B`, that induces the tableaux order.

# Examples
```jldoctest
julia> point_ordering(BracketAlgebra(6,3,[2,1,3,4,5,6]))
6-element Vector{Int64}:
 2
 1
 3
 4
 5
 6
```
"""
point_ordering(B::BracketAlgebra) = B.point_ordering

"""
    point_ordering!(B::BracketAlgebra, point_ordering::AbstractVector{<:Integer}=collect(1:B.n))

Update the point ordering of the bracket algebra `B` to `point_ordering` and return the updated bracket algebra.

Warning: changing the point ordering will not change the representation of already generated elements of `B`, which leads to wrong results when working with these elements. It is recommended to create a new bracket algebra with the desired point ordering instead. 

# Examples
```jldoctest
julia> point_ordering!(BracketAlgebra(6,3), [6,5,4,3,2,1])
Bracket algebra over P^3 on 6 points with point ordering 6 < 5 < 4 < 3 < 2 < 1 and coefficient ring Integers.
```
"""
function point_ordering!(B::BracketAlgebra, point_ordering::AbstractVector{<:Integer}=collect(1:B.n))
    B.point_ordering = point_ordering
    brackets = sort(sort.(collect(combinations(1:n(B), d(B) + 1)), lt=_lt(point_ordering)), lt=_lt(point_ordering))

    vars = AbstractAlgebra.variable_names("x#" => brackets)
    B.R, x = polynomial_ring(base_ring(B.R), vars; internal_ordering=:degrevlex)
    B.variables = Bijection{Vector{Int},typeof(x[1])}()

    for (i, bracket) in enumerate(brackets)
        B.variables[bracket] = x[i]
    end

    return B
end

# aux function that returns the less than function from the point_ordering(B)
_lt(B::BracketAlgebra) = _lt(point_ordering(B))

function Base.show(io::IO, B::BracketAlgebra)
    ordering_str = string(point_labels(B)[point_ordering(B)[1]]) * ""
    for i in 2:n(B)
        ordering_str *= " < $(point_labels(B)[point_ordering(B)[i]])"
    end
    print(io, "Bracket algebra over P^$(d(B)) on $(n(B)) points with point ordering $(ordering_str) and coefficient ring $(base_ring(B)).")
end

"""
    BracketAlgebraElem{T<:Union{RingElem,Number}} <: AbstractBracketAlgebraElem{T}

Element type for [`BracketAlgebra`](@ref) with coefficients of type `T`.
    
# Fields
- `parent::BracketAlgebra{T}`: parent object of the bracket algebra element.
- `polynomial::MPolyRingElem{T}`: polynomial in `parent.R` that represents the element.
"""
mutable struct BracketAlgebraElem{T<:Union{RingElem,Number}} <: AbstractBracketAlgebraElem{T}
    parent::BracketAlgebra{T}
    polynomial::MPolyRingElem{T}
end

function Base.show(io::IO, b::BracketAlgebraElem)
    if iszero(b)
        print(io, "0")
        return
    end

    if isone(b)
        print(io, "1")
        return
    end

    B = parent(b)

    exponents = exponent_vectors(b)
    str = ""

    for (i, exp) in enumerate(exponents)
        coefficient = coeff(b, i)

        if sign(coefficient) == -1
            str = str * " - "
        elseif i > 1
            str = str * " + "
        end

        if !(coefficient in [one(base_ring(b)), -one(base_ring(b))])
            str = str * "$(abs(coefficient))"
        end

        for (j, val) in enumerate(exp)
            if val == 0
                continue
            end
            # Retrieve the original bracket
            bracket = B.variables(gens(B.R)[j]) # Bracket is now a vector of point indices.
            # Map point indices to point labels.
            label_str = string(point_labels(B)[bracket])
            if val == 1
                str *= label_str
            else
                str *= label_str * "^$val"
            end
        end
    end

    print(io, str)
end

include("ring_interface.jl")
include("straightening.jl")

end