module BracketAlgebras

import Nemo
using Bijections
using Combinatorics

# Write your package code here.
abstract type AbstractBracketAlgebra{T} <: Nemo.Ring end
abstract type AbstractBracketAlgebraElem{T} <: Nemo.RingElem end

export AbstractBracketAlgebra, BracketAlgebra, sizyges, Tabloid, is_standard, bracket_monomial, reduced_groebner_basis!, AbstractBracketAlgebraElem, BracketAlgebraElem, point_ordering!, straighten


mutable struct BracketAlgebra{T<:Nemo.RingElement} <: AbstractBracketAlgebra{T}
    d::Int # dimension of the underlying projective space $\mathbb{P}^d$
    n::Int # number of points considered
    R::Nemo.MPolyRing{T}
    point_ordering::Vector{Int} # ordering of the n points, that induces the tableaux order. point_ordering[1] < point_ordering[2] < ... < point_ordering[n]. See Sturmfels 2008 p.81
    variables::Bijection{Vector{Int},<:Nemo.MPolyRingElem{T}}
    # ordering::DegRevLex{<:Nemo.MPolyRingElem{T}}
    # groebner_basis::Union{Nothing,Vector{<:Nemo.MPolyRingElem{T}}}

    function BracketAlgebra(n, d, point_ordering=collect(1:n), T::Type=Nemo.ZZRingElem)
        if sort(point_ordering) != collect(1:n)
            error("ordering needs to order the literals 1:$n, but got $point_ordering")
        end

        S = Nemo.parent_type(T)
        brackets = reverse!(sort(sort.(collect(combinations(1:n, d + 1)), lt=_lt(point_ordering)), lt=_lt(point_ordering))) # brackets are in the form [a,b,c,...] with a > b > c according to point_ordering. They are sorted in decreasing order wrt the tableaux order extrapolated from literal_ordering.
        vars = Nemo.AbstractAlgebra.variable_names("x#" => brackets)

        # the monomial ordering is extrapolated from vars[1] > vars[2] > vars[3]... 
        # When storing monomials the constituents are stored from biggest to smallest. So vars[2] * vars[1] is stored as vars[1]*vars[2]
        # thus deglex is the usual degrevlex (the monomial is larger, which contains the largest constituent at a larger degree), as exponent vectors are compared from left to right
        R, x = Nemo.polynomial_ring(S(), vars; internal_ordering=:deglex)

        variables = Bijection{Vector{Int},typeof(x[1])}()

        for (i, bracket) in enumerate(brackets)
            variables[bracket] = x[i]
        end

        # ordering = Groebner.DegRevLex(x)

        return new{Nemo.elem_type(S)}(d, n, R, point_ordering, variables)
    end
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

    vars = Nemo.AbstractAlgebra.variable_names("x#" => brackets)
    B.R, x = Nemo.polynomial_ring(Nemo.base_ring(B.R), vars; internal_ordering=:degrevlex)
    B.variables = Bijection{Vector{Int},typeof(x[1])}()

    for (i, bracket) in enumerate(brackets)
        B.variables[bracket] = x[i]
    end

    return B
end

function Base.show(io::IO, B::BracketAlgebra{T}) where {T}
    print(io, "Bracket algebra over P^$(d(B)) on $(n(B)) points.")
end

mutable struct BracketAlgebraElem{T<:Union{Nemo.RingElem,Number}} <: AbstractBracketAlgebraElem{T}
    parent::BracketAlgebra{T}
    polynomial::Nemo.MPolyRingElem{T}
end

function Base.show(io::IO, b::BracketAlgebraElem)
    exponents = Nemo.exponent_vectors(b)
    str = ""

    for (i, exp) in enumerate(exponents)
        coeff = Nemo.coeff(b, i)

        if Nemo.sign(coeff) == -1
            str = str * " - "
        elseif i > 1
            str = str * " + "
        end

        if !(coeff in [Nemo.one(Nemo.base_ring(b)), -Nemo.one(Nemo.base_ring(b))])
            str = str * "$coeff"
        end

        for (j, val) in enumerate(exp)
            if val == 0
                continue
            elseif val == 1
                str = str * "$(parent(b).variables(Nemo.gens(parent(b).R)[j]))"
            else
                str = str * "$(parent(b).variables(Nemo.gens(parent(b).R)[j]))" * "^$val"
            end
        end
    end

    print(io, str)
end

include("ring_interface.jl")
include("basics.jl")

end
