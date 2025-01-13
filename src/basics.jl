# aux function for the lexicographic order of vectors extended from literal_ordering
# literal_ordering[1] < literal_ordering[2] < ...
# See Sturmfels 2008 p.81
function _lt(literal_ordering::Vector{Int})
    function lt(a::Integer, b::Integer)
        indexin(a, literal_ordering)[1] < indexin(b, literal_ordering)[1]
    end
    function lt(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
        return indexin(sort(a, lt=lt), literal_ordering) < indexin(sort(b, lt=lt), literal_ordering)
    end
    return lt
end

_lt(B::BracketAlgebra) = _lt(B.literal_ordering)

function ordering!(B::BracketAlgebra, literal_ordering::AbstractVector{<:Integer}=collect(1:B.n))
    B.literal_ordering = literal_ordering
    brackets = sort(sort.(collect(combinations(1:B.n, B.d + 1)), lt=_lt(literal_ordering)), lt=_lt(literal_ordering))

    vars = Nemo.AbstractAlgebra.variable_names("x#" => brackets)
    B.R, x = Nemo.polynomial_ring(Nemo.base_ring(B.R), vars; internal_ordering=:degrevlex)
    B.variables = Bijection{Vector{Int},typeof(x[1])}()

    for (i, bracket) in enumerate(brackets)
        B.variables[bracket] = x[i]
    end

    return B
end



Base.parent(b::BracketAlgebraElem) = b.parent
Nemo.elem_type(::BracketAlgebra) = BracketAlgebraElem
Nemo.parent_type(::BracketAlgebraElem) = BracketAlgebra
Nemo.base_ring(b::BracketAlgebraElem) = Nemo.base_ring(Base.parent(b).R)
Nemo.base_ring(B::BracketAlgebra) = Nemo.base_ring(B.R)

Base.one(B::BracketAlgebra) = BracketAlgebraElem(B, one(B.R))
Base.zero(B::BracketAlgebra) = BracketAlgebraElem(B, zero(B.R))
Base.one(b::BracketAlgebraElem) = one(parent(b))
Base.zero(b::BracketAlgebraElem) = zero(parent(b))


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

import Base.==
function ==(a::BracketAlgebraElem, b::BracketAlgebraElem)
    return (parent(a).n == parent(b).n) && (parent(a).d == parent(b).d) && (parent(a).literal_ordering == parent(b).literal_ordering) && straighten(a - b).polynomial == zero(parent(b).R)
end


function Nemo.evaluate(b::BracketAlgebraElem{T}, A::Vector{T}) where {T<:Union{Nemo.RingElem,Number}}
    Nemo.evaluate(b.polynomial, A)
end

function Nemo.evaluate(b::BracketAlgebraElem{T}, A::Vector{U}) where {T<:Union{Nemo.RingElem,Number},U<:Integer}
    Nemo.evaluate(b.polynomial, A)
end

function Nemo.evaluate(b::BracketAlgebraElem{T}, coordinization::AbstractMatrix{<:Union{Nemo.RingElem,Number}}) where {T<:Union{Nemo.RingElem,Number}}
    bracks = brackets(b)
    A = map(x -> x in bracks ? Nemo.det(T.(hcat(transpose(coordinization[:, x]), ones(parent(b).d + 1, 1)))) : 0, sort(collect(keys(parent(b).variables)), rev=true))
    Nemo.evaluate(b.polynomial, A)
end

mutable struct Tabloid
    matrix::Matrix{Int}
    ordering::Vector{Int}

    """
    Tabloid(matrix::AbstractMatrix{<:Integer}, ordering::AbstractVector{<:Integer}=collect(1:maximum(matrix)))

    Calculate the tabloid whose rows correspond to the rows of matrix with the ordering.
        Example:
        matrix = [4 2 1; 3 1 4], ordering = [3,1,2,4]
        => result: [3 1 4; 1 2 4] 
    """
    function Tabloid(matrix::AbstractMatrix{<:Integer}, ordering::AbstractVector{<:Integer}=collect(1:maximum(matrix)))
        rows = [matrix[i, :] for i in 1:size(matrix)[1]]
        sort!.(rows, by=(x -> indexin(x, ordering)[1]))
        sort!(rows, by=(row -> indexin(row, ordering)))
        return new(transpose(hcat(rows...)), ordering)
    end
end

Base.show(io::IO, t::Tabloid) = show(io, t.matrix)

function Tabloid(rows::AbstractVector{<:AbstractVector{<:Integer}}, ordering::AbstractVector{<:Integer}=collect(1:maximum(vcat(rows...))))
    return Tabloid(transpose(hcat(rows...)), ordering)
end

Matrix(t::Tabloid) = t.matrix

Base.vcat(t1::Tabloid, t2::Tabloid) = t1.ordering == t2.ordering ? Tabloid(t1.ordering, vcat(t1.matrix, t2.matrix)) : error("Tabloids need to have same number of columns.")


"""
    tabloid(b::BracketAlgebraElem, ordering::Vector{Int}=collect(1:parent(b).n))

Return the tabloid corresponding to the bracket monomial b ordered by ordering. 
ordering determines the order in which 

Example:
b = [1,2,3][3,4,5], ordering = [4,3,1,2,5]
result: [4 3 5; 3 1 2]
"""
function Tabloid(b::BracketAlgebraElem)
    # see Sturmfels 2008, page 81 on how to build the tabloids
    if length(b) > 1
        error("Only tabloids of bracket monomials can be calculated.")
    end

    exps = collect(Nemo.exponent_vectors(b))[1]
    # all brackets that appear as rows 
    rows = [(repeat(parent(b).variables(Nemo.gens(parent(b).R)[i]), exps[i]) for i in eachindex(exps))...]
    filter!(row -> length(row) > 0, rows)
    return Tabloid(rows, parent(b).literal_ordering)
end

"""
    standard_violation(t::Tabloid)

Return the index of the first violation to the standardness of t, i.e. the first index where t[i,j] > t[i+1, j] with regard to the ordering. Otherwise return nothing.
Example:
t = [1 2 3; 1 4 5; 1 5 6; 2 3 4] with ordering [1,2,3,4,5,6] => return (3,2)
t = [1 2 3; 1 2 4] with ordering [1,2,3,4] => return nothing
"""
function standard_violation(t::Tabloid)
    if size(t.matrix)[1] == 1
        return nothing
    end

    return findfirst([_lt(t.ordering)(t.matrix[row+1, col], t.matrix[row, col]) for row in 1:size(t.matrix)[1]-1, col in 1:size(t.matrix)[2]])
end

standard_violation(b::BracketAlgebraElem, ordering::Vector{<:Integer}=collect(1:parent(b).n)) = standard_violation(Tabloid(b, ordering))

is_standard(t::Tabloid) = isnothing(standard_violation(t))
is_standard(b::BracketAlgebraElem) = is_standard(Tabloid(b))

function straightening_sizyge(α::Vector{<:Integer}, β::Vector{<:Integer}, γ::Vector{<:Integer}, B::BracketAlgebra)
    s = length(α) + 1
    d = B.d
    @assert length(β) == d + 2 "β needs to have length d + 2, but got $(length(β))."
    @assert length(γ) == d + 1 - s "γ needs to have length d + 1 - s, but got $(length(γ))."

    return sum(Nemo.sign(Nemo.Perm(vcat(setdiff(collect(1:d+2), τ), τ))) * B(vcat(α, β[setdiff(collect(1:d+2), τ)])) * B(vcat(β[τ], γ)) for τ in combinations(1:d+2, s))
end

"""
    straighten(b::BracketAlgebraElem, ordering::AbstractVector{<:Integer} = collect(1:parent(b).n); chek::Bool = true)

Perform the straightening algorithm to the BracketAlgebraElem b with monomial ordering induced by ordering[1] < ordering[2] < ...

For details see Sturmfels 2008, chapter 3.1, Handbook of Geometric Constraint System Principles, chapater 4.3
"""
function straighten(b::BracketAlgebraElem)
    first_nonstandard_ind = findfirst(mon -> !is_standard(mon), collect(Nemo.monomials(b)))
    if isnothing(first_nonstandard_ind)
        return b
    end

    first_nonstandard = collect(Nemo.monomials(b))[first_nonstandard_ind]
    t = Tabloid(first_nonstandard)
    (r, s) = Tuple(standard_violation(t))

    mat = Matrix(t)
    α = mat[r, 1:s-1]
    β = vcat(mat[r, s:end], mat[r+1, 1:s])
    γ = mat[r+1, s+1:end]

    sizyge = straightening_sizyge(α, β, γ, parent(b))
    return straighten(b - (Nemo.coeff(b, first_nonstandard_ind)) * prod(parent(b)(mat[i, :]) for i in setdiff(1:size(mat)[1], [r, r + 1]); init=one(parent(b))) * sizyge)
end

# function to find atomic extensors of BracketAlgebraElem as in Sturmfels 2008, Alg. 3.5.6 step 1
"""
    atomic_extensors(b::BracketAlgebraElem)

Find all atomic extensors of the bracket algebra elem b. These are the equivalence classes of the relation p1 ~ p2 iff substituting p1 = p2 in b results in the zero element of the bracket algebra.
"""
function atomic_extensors(b::BracketAlgebraElem)
    B = parent(b)

    extensors = Vector{Int}[]

    exps = Nemo.exponent_vectors(b)
    coeffs = Nemo.coefficients(b)
    function ~(p, q)

    end
end