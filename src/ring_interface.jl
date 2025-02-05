

Base.parent(b::BracketAlgebraElem) = b.parent
Nemo.elem_type(::BracketAlgebra) = BracketAlgebraElem
Nemo.parent_type(::BracketAlgebraElem) = BracketAlgebra
Nemo.base_ring(b::BracketAlgebraElem) = Nemo.base_ring(Base.parent(b).R)
Nemo.base_ring(B::BracketAlgebra) = Nemo.base_ring(B.R)

Base.one(B::BracketAlgebra) = BracketAlgebraElem(B, one(B.R))
Base.zero(B::BracketAlgebra) = BracketAlgebraElem(B, zero(B.R))
Base.one(b::BracketAlgebraElem) = one(parent(b))
Base.zero(b::BracketAlgebraElem) = zero(parent(b))

############################################################################################################
# Constructors for bracket algebra elements
############################################################################################################
(B::BracketAlgebra)(A::Vector{T}, m::Vector{Vector{Int}}) where {T<:Nemo.RingElem} = BracketAlgebraElem(B, B.R(A, m))
(B::BracketAlgebra)(p::Nemo.MPolyRingElem) = BracketAlgebraElem(B, p)

# Bracket expression from array: B([1,2,3,4]) = [1,2,3,4]
(B::BracketAlgebra)(bracket::Vector{<:Integer}) = length(unique(bracket)) == length(bracket) ? Nemo.sign(Nemo.Perm(Int.(indexin(bracket, sort(bracket, lt=_lt(B))))))BracketAlgebraElem(B, B.variables[sort(bracket, lt=_lt(B))]) : zero(B)

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

    exponent_vecs_A = collect(Nemo.exponent_vectors(a.polynomial))

    # determine the vector of indices of the variables in Nemo.gens(A.R) as variables in Nemo.gens(B.R) and a vector of the sign changes when translating the variables
    translate = zeros(Int, length(Nemo.gens(A.R)))
    sign_changes = ones(Int, length(Nemo.gens(A.R)))
    for (i, var) in enumerate(Nemo.gens(A.R))
        bracket_a = A.variables(var) # vector of ints corresponding to var
        b = B(bracket_a) # bracket algebra element of B corresponding to var (is a bracket algebra monomial of degree 1)

        sign_changes[i] = Nemo.coeff(b.polynomial, 1) # sign change is the sign of b

        var_B = (Nemo.gens(B.R)[Nemo.exponent_vector(b.polynomial, 1).>0])[1] # variable in B corresponding to var
        translate[i] = indexin([var_B], Nemo.gens(B.R))[1] # index of var_B in Nemo.gens(B.R)
    end

    exponent_vecs_B = [zeros(Int, length(Nemo.gens(B.R))) for _ in 1:length(exponent_vecs_A)]

    for (i, exp) in enumerate(exponent_vecs_B)
        exp[translate] = exponent_vecs_A[i]
    end

    coeffs_B = [prod(sign_changes .^ exp) * Nemo.coeff(a.polynomial, i) for (i, exp) in enumerate(exponent_vecs_A)]

    return B(coeffs_B, exponent_vecs_B)
end

Base.length(b::BracketAlgebraElem) = Nemo.length(b.polynomial)
Nemo.degrees(b::BracketAlgebraElem) = Nemo.degrees(b.polynomial)
Nemo.total_degree(b::BracketAlgebraElem) = Nemo.total_degree(b.polynomial)
Nemo.coefficients(b::BracketAlgebraElem) = Nemo.coefficients(b.polynomial)
Nemo.monomials(b::BracketAlgebraElem) = (parent(b)(p) for p in Nemo.monomials(b.polynomial))
Nemo.terms(b::BracketAlgebraElem) = (parent(b)(p) for p in Nemo.terms(b.polynomial))
Nemo.exponent_vectors(b::BracketAlgebraElem) = Nemo.exponent_vectors(b.polynomial)
Nemo.coeff(b::BracketAlgebraElem, n::Int) = Nemo.coeff(b.polynomial, n)
Nemo.coeff(b::BracketAlgebraElem, exps::Vector{Int}) = Nemo.coeff(b.polynomial, exps)
Nemo.monomial(b::BracketAlgebraElem, n::Int) = parent(b)(Nemo.monomial(b.polynomial, n))
Nemo.exponent_vector(b::BracketAlgebra, n::Int) = Nemo.exponent_vector(b.polynomial, n)
Nemo.term(b::BracketAlgebraElem, n::Int) = parent(b)(Nemo.term(b.polynomial, n))
Nemo.leading_term(b::BracketAlgebraElem) = parent(b)(Nemo.leading_term(b.polynomial))
Nemo.leading_monomial(b::BracketAlgebraElem) = parent(b)(Nemo.leading_monomial(b.polynomial))
Nemo.leading_coefficient(b::BracketAlgebraElem) = Nemo.leading_coefficient(b.polynomial)

Nemo.factor(b::BracketAlgebraElem) = Nemo.factor(b.polynomial)

# return all brackets that appear in b as arrays
brackets(b::BracketAlgebraElem) = sort(collect(keys(parent(b).variables)), rev=true)[sum(Nemo.exponent_vectors(b)).>0]

Base.:*(n::Integer, b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, n * b.polynomial)
Base.:*(c::T, b::BracketAlgebraElem{T}) where {T<:Nemo.RingElem} = parent(b)(c * b.polynomial)
Base.:*(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial * b.polynomial)
Base.:+(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial + b.polynomial)
Base.:-(a::BracketAlgebraElem, b::BracketAlgebraElem) = BracketAlgebraElem(a.parent, a.polynomial - b.polynomial)
Base.:-(b::BracketAlgebraElem) = BracketAlgebraElem(b.parent, -b.polynomial)
Base.:^(b::BracketAlgebraElem, n::Int) = BracketAlgebraElem(b.parent, b.polynomial^n)
Base.:>(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial > b.polynomial
Base.:<(a::BracketAlgebraElem, b::BracketAlgebraElem) = a.polynomial < b.polynomial