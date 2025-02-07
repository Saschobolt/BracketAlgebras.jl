############################################################################################################
# AbstractAlgebra Ring Interface. See https://nemocas.github.io/AbstractAlgebra.jl/stable/ring_interface/
############################################################################################################

############################################################################################################
# Data type and parent object methods
############################################################################################################
parent(b::BracketAlgebraElem) = b.parent
elem_type(::BracketAlgebra) = BracketAlgebraElem
parent_type(::BracketAlgebraElem) = BracketAlgebra
base_ring(b::BracketAlgebraElem) = base_ring(parent(b).R)
base_ring(B::BracketAlgebra) = base_ring(B.R)
is_domain_type(::Type{BracketAlgebraElem}) = false
is_exact_type(::Type{BracketAlgebraElem{T}}) where {T} = is_exact_type(T)
Base.deepcopy_internal(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, deepcopy(b.polynomial))

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
# Constructors for bracket algebra elements
############################################################################################################
(B::BracketAlgebra)(A::Vector{T}, m::Vector{Vector{Int}}) where {T<:RingElem} = BracketAlgebraElem(B, B.R(A, m))
(B::BracketAlgebra)(p::MPolyRingElem) = BracketAlgebraElem(B, p)
(B::BracketAlgebra)(n::Integer) = BracketAlgebraElem(B, B.R(n))
(B::BracketAlgebra{T})(n::T) where {T<:RingElem} = BracketAlgebraElem(B, B.R(n))

# Bracket expression from array: B([1,2,3,4]) = [1,2,3,4]
(B::BracketAlgebra)(bracket::Vector{<:Integer}) = length(unique(bracket)) == length(bracket) ? sign(Perm(Int.(indexin(bracket, sort(bracket, lt=_lt(B))))))BracketAlgebraElem(B, B.variables[sort(bracket, lt=_lt(B))]) : zero(B)

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