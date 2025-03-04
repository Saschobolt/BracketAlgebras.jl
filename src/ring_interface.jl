############################################################################################################
# AbstractAlgebra Ring Interface. See https://nemocas.github.io/AbstractAlgebra.jl/stable/ring_interface/
############################################################################################################
import AbstractAlgebra: parent_type, elem_type, base_ring, base_ring_type, parent, is_domain_type, is_exact_type, canonical_unit, isequal, divexact, zero!, mul!, add!, get_cached!, is_unit, characteristic,
    Ring, RingElem, RingElement, MPolyRing, MPolyRingElem,
    degrees, total_degree, coefficients, monomials, terms, exponent_vectors, coeff, monomial, exponent_vector, term, leading_term, leading_monomial, leading_coefficient, factor, evaulate

import Base: rand, one, zero, length, iszero, isone, promote_rule, deepcopy, *, +, -, ^, >, <, ==, hash, deepcopy_internal
############################################################################################################
# Data type and parent object methods
############################################################################################################
parent_type(::Type{BracketAlgebraElem{T}}) where {T<:RingElement} = BracketAlgebra{T}

elem_type(::Type{BracketAlgebra{T}}) where {T<:RingElement} = BracketAlgebraElem{T}

base_ring_type(::Type{BracketAlgebra{T}}) where {T<:RingElement} = parent_type(T)

base_ring(b::BracketAlgebraElem) = base_ring(parent(b).R)

base_ring(B::BracketAlgebra) = base_ring(B.R)

parent(b::BracketAlgebraElem) = b.parent

is_domain_type(::Type{BracketAlgebraElem}) = false

is_exact_type(::Type{BracketAlgebraElem{T}}) where {T} = is_exact_type(T)

deepcopy_internal(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, deepcopy(b.polynomial))

function hash(b::BracketAlgebraElem, h::UInt)
    r = 0xc2fa78ade36f978d
    return xor(r, hash(b.polynomial, h))
end

############################################################################################################
# Constructors for bracket algebra elements
############################################################################################################
(B::BracketAlgebra)() = zero(B)

(B::BracketAlgebra)(a::Integer) = BracketAlgebraElem(B, B.R(a))

(B::BracketAlgebra{T})(a::T) where {T<:RingElement} = BracketAlgebraElem(B, B.R(a))

(B::BracketAlgebra{T})(A::Vector{T}, m::Vector{<:Vector{<:Integer}}) where {T<:RingElement} = BracketAlgebraElem(B, B.R(A, m))

(B::BracketAlgebra)(p::MPolyRingElem) = BracketAlgebraElem(B, p)

(B::BracketAlgebra{T})(n::T) where {T<:RingElem} = BracketAlgebraElem(B, B.R(n))

# Bracket expression from array of point indices: B([1,2,3,4]) = [1,2,3,4]
(B::BracketAlgebra)(bracket::Vector{<:Integer}) = length(unique(bracket)) == length(bracket) ? sign(Perm(Int.(indexin(bracket, sort(bracket, lt=_lt(B))))))BracketAlgebraElem(B, B.variables[sort(bracket, lt=_lt(B))]) : zero(B)

# Constructor for bracket algebra elements by point labels. E.g. B = BracketAlgebra(4,2, point_labels = ["a", "b", "c", "d"]) => B(["b", "a","c"]) == B([2,1,3]) == -["a", "b", "c"]
function (B::BracketAlgebra)(bracket::AbstractVector)
    if !issubset(bracket, point_labels(B))
        error("`bracket` ($bracket) doesn't match point labels of bracket algebra.")
    end

    return B(Vector{Int}(indexin(bracket, point_labels(B))))
end

# Bracket polynomial from array of array of arrays. They encode the bracket polynomial as a sum of monomials. B([[[1,2], [3,4]], [2,3]]) = [1,2]*[3,4] + [2,3]
(B::BracketAlgebra)(A::Vector{<:Vector{<:Vector{<:Integer}}}) = sum(prod(B(bracket) for bracket in monomial) for monomial in A)

# rewrite the element a of bracket algebra A as an element of the bracket algebra B 
function (B::BracketAlgebra)(a::BracketAlgebraElem)
    A = parent(a)
    if A.d != B.d
        error("parent of a and B have to model same projective space, but B.d = $(B.d) and parent(a).d = $(A.d).")
    end
    if A.n > B.n
        error("parent of a has to model less points than B, but B.n = $(B.n) and parent(a).n = $(A.n).")
    end

    exponent_vecs_A = collect(exponent_vectors(a.polynomial))

    # determine the vector of indices of the variables in gens(A.R) as variables in gens(B.R) and a vector of the sign changes when translating the variables
    translate = zeros(Int, length(gens(A.R)))
    sign_changes = ones(Int, length(gens(A.R)))
    for (i, var) in enumerate(gens(A.R))
        bracket_a = A.variables(var) # vector of ints corresponding to var
        b = B(bracket_a) # bracket algebra element of B corresponding to var (is a bracket algebra monomial of degree 1)

        sign_changes[i] = coeff(b.polynomial, 1) # sign change is the sign of b

        var_B = (gens(B.R)[exponent_vector(b.polynomial, 1).>0])[1] # variable in B corresponding to var
        translate[i] = indexin([var_B], gens(B.R))[1] # index of var_B in gens(B.R)
    end

    exponent_vecs_B = [zeros(Int, length(gens(B.R))) for _ in 1:length(exponent_vecs_A)]

    for (i, exp) in enumerate(exponent_vecs_B)
        exp[translate] = exponent_vecs_A[i]
    end

    coeffs_B = [prod(sign_changes .^ exp) * coeff(a.polynomial, i) for (i, exp) in enumerate(exponent_vecs_A)]

    return B(coeffs_B, exponent_vecs_B)
end

############################################################################################################
# Basic Manipulation of rings and elements
############################################################################################################
one(B::BracketAlgebra) = BracketAlgebraElem(B, one(B.R))

zero(B::BracketAlgebra) = BracketAlgebraElem(B, zero(B.R))

one(b::BracketAlgebraElem) = one(parent(b))

zero(b::BracketAlgebraElem) = zero(parent(b))

iszero(b::BracketAlgebraElem) = iszero(b.polynomial)

isone(b::BracketAlgebraElem) = isone(b.polynomial)

canonical_unit(b::BracketAlgebraElem) = parent(b)(canonical_unit(b.polynomial))

############################################################################################################
# Promotion rules
############################################################################################################
promote_rule(::Type{BracketAlgebraElem{T}}, ::Type{BracketAlgebraElem{T}}) where {T<:RingElement} = BracketAlgebraElem{T}

function promote_rule(::Type{BracketAlgebraElem{T}}, ::Type{U}) where {T<:RingElement,U<:RingElement}
    promote_rule(T, U) == T ? BracketAlegbraElem{T} : Union{}
end

############################################################################################################
# functions for BracketAlgebraElem as polynomial
############################################################################################################
length(b::BracketAlgebraElem) = length(b.polynomial)

degrees(b::BracketAlgebraElem) = degrees(b.polynomial)

total_degree(b::BracketAlgebraElem) = total_degree(b.polynomial)

coefficients(b::BracketAlgebraElem) = coefficients(b.polynomial)

monomials(b::BracketAlgebraElem) = (parent(b)(p) for p in monomials(b.polynomial))

terms(b::BracketAlgebraElem) = (parent(b)(p) for p in terms(b.polynomial))

exponent_vectors(b::BracketAlgebraElem) = exponent_vectors(b.polynomial)

coeff(b::BracketAlgebraElem, n::Int) = coeff(b.polynomial, n)

coeff(b::BracketAlgebraElem, exps::Vector{Int}) = coeff(b.polynomial, exps)

monomial(b::BracketAlgebraElem, n::Int) = parent(b)(monomial(b.polynomial, n))

exponent_vector(b::BracketAlgebra, n::Int) = exponent_vector(b.polynomial, n)

term(b::BracketAlgebraElem, n::Int) = parent(b)(term(b.polynomial, n))

leading_term(b::BracketAlgebraElem) = parent(b)(leading_term(b.polynomial))

leading_monomial(b::BracketAlgebraElem) = parent(b)(leading_monomial(b.polynomial))

leading_coefficient(b::BracketAlgebraElem) = leading_coefficient(b.polynomial)

factor(b::BracketAlgebraElem) = factor(b.polynomial)

# return all brackets that appear in b as arrays
brackets(b::BracketAlgebraElem) = sort(collect(keys(parent(b).variables)), rev=true)[sum(exponent_vectors(b)).>0]

function evaluate(b::BracketAlgebraElem{T}, A::Vector{T}) where {T<:Union{RingElem,Number}}
    evaluate(b.polynomial, A)
end

function evaluate(b::BracketAlgebraElem{T}, A::Vector{U}) where {T<:Union{RingElem,Number},U<:Integer}
    evaluate(b.polynomial, A)
end

function evaluate(b::BracketAlgebraElem{T}, coordinization::AbstractMatrix{<:Union{RingElem,Number}}) where {T<:Union{RingElem,Number}}
    bracks = brackets(b)
    A = map(x -> x in bracks ? det(T.(hcat(transpose(coordinization[:, x]), ones(parent(b).d + 1, 1)))) : 0, sort(collect(keys(parent(b).variables)), rev=true))
    evaluate(b.polynomial, A)
end

############################################################################################################
# Arithmetics and comparisons (no ==, this is provided in straightening.jl)
############################################################################################################
Base.:*(n::Integer, b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, n * b.polynomial)

Base.:*(c::T, b::BracketAlgebraElem{T}) where {T<:RingElem} = parent(b)(c * b.polynomial)

Base.:*(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial * b.polynomial)

Base.:+(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial + b.polynomial)

Base.:-(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial - b.polynomial)

Base.:-(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, -b.polynomial)

Base.:^(b::BracketAlgebraElem, n::Int) = BracketAlgebraElem(b.parent, b.polynomial^n)

Base.:>(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial > b.polynomial

Base.:<(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial < b.polynomial

function ==(a::BracketAlgebraElem, b::BracketAlgebraElem)
    return (parent(a).n == parent(b).n) && (parent(a).d == parent(b).d) && (point_ordering(parent(a)) == point_ordering(parent(b))) && straighten(a - b).polynomial == zero(parent(b).R)
end

############################################################################################################
# Random generation
############################################################################################################

# TODO

