module BracketAlgebras

using AbstractAlgebra
# import Base: show, +, -, *, ^, ==, <, >, inv, isone, iszero, one, zero, rand, deepcopy_internal, hash
using Bijections
using Combinatorics

abstract type AbstractBracketAlgebra{T} <: Ring end
abstract type AbstractBracketAlgebraElem{T} <: RingElem end

export AbstractBracketAlgebra, BracketAlgebra, sizyges, Tabloid, is_standard, bracket_monomial, reduced_groebner_basis!, AbstractBracketAlgebraElem, BracketAlgebraElem, point_ordering!, straighten

abstract type AbstractBracketAlgebra{T} <: Ring end
abstract type AbstractBracketAlgebraElem{T} <: RingElem end

export AbstractBracketAlgebra, BracketAlgebra, sizyges, Tabloid, is_standard, bracket_monomial, reduced_groebner_basis!, AbstractBracketAlgebraElem, BracketAlgebrasElem, point_ordering, point_ordering!, point_labels, point_labels!

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
    point_ordering::Vector{Int} # ordering of the n points, that induces the tableaux order. point_ordering[1] < point_ordering[2] < ... < point_ordering[n]. See Sturmfels 2008 p.81
    variables::Bijection{Vector{Int},<:MPolyRingElem{T}}
    point_labels::Vector  # New: point labels for each point


    """
        BracketAlgebra(n, d, point_ordering=collect(1:n), T::Type=Int; point_labels::Vector=collect(1:n))

    Construct a bracket algebra over the projective space P^d on n points with the given point ordering and point labels.

    # Examples
    ```jldoctest
    julia> B = BracketAlgebra(4, 2, [1,2,3,4], point_labels = ["a","b","c","d"])
    Bracket algebra over P^2 on 4 points with point ordering a < b < c < d.
    ```

    ```jldoctest
    B = BracketAlgebra(4, 2, [2,1,3,4], point_labels = ["a","b","c","d"])
    Bracket algebra over P^2 on 4 points with point ordering b < a < c < d.
    ```
    """
    function BracketAlgebra(n, d, point_ordering=collect(1:n), T::Type=Int; point_labels::Vector=collect(1:n))
        if sort(point_ordering) != collect(1:n)
            error("ordering needs to order the points 1:$n, but got $point_ordering")
        end
        if length(point_labels) != n
            error("point_labels must be of length $n")
        end

        S = parent_type(T)
        brackets = sort(sort.(collect(combinations(1:n, d + 1)), lt=_lt(point_ordering)), lt=_lt(point_ordering)) # brackets are in the form [a,b,c,...] with a < b < c according to point_ordering. They are sorted in decreasing order wrt the tableaux order extrapolated from point_ordering.
        vars = AbstractAlgebra.variable_names("x#" => brackets)

        # the monomial ordering is extrapolated from vars[1] > vars[2] > vars[3]... As the brackets are sorted in increasing order wrt the tableaux order, the variables need to be reversed.
        R, x = polynomial_ring(S(), reverse(vars), internal_ordering=:degrevlex)
        reverse!(x) # now x is in the same order as brackets

        variables = Bijection{Vector{Int},typeof(x[1])}()

        for (i, bracket) in enumerate(brackets)
            variables[bracket] = x[i]
        end

        # ordering = Groebner.DegRevLex(x)

        return new{elem_type(S)}(d, n, R, point_ordering, variables, point_labels)
    end
end

# New accessor for point_labels.
point_labels(B::BracketAlgebra) = B.point_labels

# New modifier for point_labels.
function point_labels!(B::BracketAlgebra, new_labels::AbstractVector)
    if length(new_labels) != B.n
        error("new_labels must be of length $(B.n)")
    end
    B.point_labels = new_labels
    return B
end

"""
    d(B::BracketAlgebra)

Return the dimension of the projective space \$\\mathbb{P}^d\$ underlying the bracket algebra B.
"""
d(B::BracketAlgebra) = B.d

"""
    n(B::BracketAlgebra)

Return the number of points considered in the bracket algebra B.
"""
n(B::BracketAlgebra) = B.n

"""
    point_ordering(B::BracketAlgebra)

Return the ordering of the n points in the bracket algebra B, that induces the tableaux order.
"""
point_ordering(B::BracketAlgebra) = B.point_ordering

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

function Base.show(io::IO, B::BracketAlgebra)
    ordering_str = string(point_labels(B)[point_ordering(B)[1]]) * ""
    for i in 2:n(B)
        ordering_str *= " < $(point_labels(B)[point_ordering(B)[i]])"
    end
    print(io, "Bracket algebra over P^$(d(B)) on $(n(B)) points with point ordering $(ordering_str).")
end

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
            str = str * "$coefficient"
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
include("basics.jl")

end