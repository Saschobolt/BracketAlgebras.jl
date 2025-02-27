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

    exps = collect(exponent_vectors(b))[1]
    # all brackets that appear as rows 
    rows = [(repeat(parent(b).variables(gens(parent(b).R)[i]), exps[i]) for i in eachindex(exps))...]
    filter!(row -> length(row) > 0, rows)
    return Tabloid(rows, point_ordering(parent(b)))
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

    return sum(sign(Perm(vcat(setdiff(collect(1:d+2), τ), τ))) * B(vcat(α, β[setdiff(collect(1:d+2), τ)])) * B(vcat(β[τ], γ)) for τ in combinations(1:d+2, s))
end

"""
    straighten(b::BracketAlgebraElem, ordering::AbstractVector{<:Integer} = collect(1:parent(b).n); chek::Bool = true)

Perform the straightening algorithm to the BracketAlgebraElem b with monomial ordering induced by ordering[1] < ordering[2] < ...

For details see Sturmfels 2008, chapter 3.1, Handbook of Geometric Constraint System Principles, chapater 4.3
"""
function straighten(b::BracketAlgebraElem)
    first_nonstandard_ind = findfirst(mon -> !is_standard(mon), collect(monomials(b)))
    if isnothing(first_nonstandard_ind)
        return b
    end

    first_nonstandard = collect(monomials(b))[first_nonstandard_ind]
    t = Tabloid(first_nonstandard)
    (r, s) = Tuple(standard_violation(t))

    mat = Matrix(t)
    α = mat[r, 1:s-1]
    β = vcat(mat[r, s:end], mat[r+1, 1:s])
    γ = mat[r+1, s+1:end]

    sizyge = straightening_sizyge(α, β, γ, parent(b))
    return straighten(b - (coeff(b, first_nonstandard_ind)) * prod(parent(b)(mat[i, :]) for i in setdiff(1:size(mat)[1], [r, r + 1]); init=one(parent(b))) * sizyge)
end

# function to find atomic extensors of BracketAlgebraElem as in Sturmfels 2008, Alg. 3.5.6 step 1
"""
    atomic_extensors(b::BracketAlgebraElem)

Find all atomic extensors of the bracket algebra elem b. These are the equivalence classes of the relation p1 ~ p2 iff substituting p1 = p2 in b results in the zero element of the bracket algebra.
"""
function atomic_extensors(b::BracketAlgebraElem)
    B = parent(b)

    extensors = Vector{Int}[]

    exps = exponent_vectors(b)
    coeffs = coefficients(b)
    function ~(p, q)

    end
end