module BracketAlgebras

import Nemo
using Bijections
using Combinatorics

# Write your package code here.
abstract type AbstractBracketAlgebra <: Nemo.Ring end
abstract type AbstractBracketAlgebraElem <: Nemo.RingElem end

export AbstractBracketAlgebra, BracketAlgebra, sizyges, Tabloid, is_standard, bracket_monomial, reduced_groebner_basis!, AbstractBracketAlgebraElem, BracketAlgebraElem, ordering!, straighten


struct BracketAlgebra{T<:Union{Nemo.RingElem,Number}} <: AbstractBracketAlgebra
    d::Int # dimension of the modelled projective space $\mathbb{P}^d$
    n::Int # number of points in the projective space
    R::Nemo.MPolyRing{T}
    literal_ordering::Vector{Int} # ordering of the n points, that induces the tableaux order. literal_ordering[1] < literal_ordering[2] < ... < literal_ordering[n]. See Sturmfels 2008 p.81
    variables::Bijection{Vector{Int},<:Nemo.MPolyRingElem{T}}
    # ordering::DegRevLex{<:Nemo.MPolyRingElem{T}}
    # groebner_basis::Union{Nothing,Vector{<:Nemo.MPolyRingElem{T}}}

    function BracketAlgebra(n, d, literal_ordering=collect(1:n), T::Type=Nemo.ZZRingElem)
        if sort(literal_ordering) != collect(1:n)
            error("ordering needs to order the literals 1:$n, but got $literal_ordering")
        end

        S = Nemo.parent_type(T)
        brackets = reverse!(sort(sort.(collect(combinations(1:n, d + 1)), lt=_lt(literal_ordering)), lt=_lt(literal_ordering))) # brackets are in the form [a,b,c,...] with a > b > c according to literal_ordering. They are sorted in decreasing order wrt the tableaux order extrapolated from literal_ordering.
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

        return new{Nemo.elem_type(S)}(d, n, R, literal_ordering, variables)
    end
end

mutable struct BracketAlgebraElem{T<:Union{Nemo.RingElem,Number}} <: AbstractBracketAlgebraElem
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

include("basics.jl")

end
