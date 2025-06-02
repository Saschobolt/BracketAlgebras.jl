"""
    Tabloid

A tabloid with ordered rows. 

# Fields
-`matrix::Matrix{Int}`: Matrix representation of the tabloid.
-`ordering::Vector{Int}`: Ordering of the entries of the tabloids. `ordering[1]` < `ordering[2]` < ...
"""
mutable struct Tabloid
    matrix::Matrix{Int}
    ordering::Vector{Int}

    @doc """
    Tabloid(matrix::AbstractMatrix{<:Integer}, ordering::AbstractVector{<:Integer}=collect(1:maximum(matrix)))

    Return the tabloid whose rows correspond to the rows of matrix with the given ordering. The entries of the rows are sorted wrt the ordering and then the rows are ordered lexicographically (tableaux order).

    # Examples
    ```jldoctest
    julia> Tabloid([4 2 1; 3 1 4], [3,1,2,4])
    [3 1 4; 1 2 4]
    ```

    ```jldoctest
    julia> Tabloid([4 2 1; 3 1 4], [3,2,4,1])
    [3 4 1; 2 4 1]
    ```
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
    Tabloid(b::BracketAlgebraElem)

Return the tabloid corresponding to the bracket monomial `b`.

# Examples
```jldoctest
julia> B = BracketAlgebra(5,2)
Bracket algebra over P^2 on 5 points with point ordering 1 < 2 < 3 < 4 < 5 and coefficient ring Integers.

julia> b = B([1,2,3]) * B([1,2,5])
[1, 2, 5][1, 2, 3]

julia> Tabloid(b)
[1 2 3; 1 2 5]
```
"""
function Tabloid(b::BracketAlgebraElem)
    # see Sturmfels 2008, page 81 on how to build the tabloids
    if length(b) > 1
        error("Only tabloids of bracket monomials can be calculated.")
    end

    exps = collect(exponent_vectors(b))[1]
    # all brackets that appear as rows 
    rows = Vector{Int}[]
    for (i, e) in enumerate(exps)
        for k in 1:e
            push!(rows, parent(b).variables(gens(parent(b).R)[i]))
        end
    end
    return Tabloid(rows, point_ordering(parent(b)))
end

"""
    standard_violation(t::Tabloid)

Return the index of the first violation to the standardness of the tabloid `t`, i.e. the first index where `t[i,j]` > `t[i+1, j]` with regard to the ordering. Otherwise return `nothing`.

# Examples: 
```jldoctest
julia> standard_violation(Tabloid([1 2 3; 1 4 5; 1 5 6; 2 3 4], [1,2,3,4,5,6]))
CartesianIndex(3, 2)
```

```jldoctest
julia> standard_violation(Tabloid([1 2 3; 1 2 4], [1,2,3,4]))

```
"""
function standard_violation(t::Tabloid)
    if size(t.matrix)[1] == 1
        return nothing
    end

    return findfirst([_lt(t.ordering)(t.matrix[row+1, col], t.matrix[row, col]) for row in 1:size(t.matrix)[1]-1, col in 1:size(t.matrix)[2]])
end
Vector
"""
    standard_violation(b::BracketAlgebraElem)

Return the index of the first standard violation of the tabloid corresponding to `b`.
"""
standard_violation(b::BracketAlgebraElem) = standard_violation(Tabloid(b))

"""
    is_standard(t::Tabloid)

Return whether the tabloid `t` is standard. See also [`standard_violation`](@ref).
"""
is_standard(t::Tabloid) = isnothing(standard_violation(t))

"""
    is_standard(b::BracketAlgebraElem)

Return whether the tabloid correspondong to the bracket algebra element `b` is standard. See also [`standard_violation`](@ref).
"""
is_standard(b::BracketAlgebraElem) = is_standard(Tabloid(b))


"""
    straightening_sizyge(α::Vector{<:Integer}, β::Vector{<:Integer}, γ::Vector{<:Integer}, B::BracketAlgebra)

Return the straightening_sizyge corresponding to `α`, `β`, `γ` as an element of the bracket algebra `B`.
The lenghts of the vectors have to fulfill `length(β)==d(B)+2` and `length(γ)==d(B)-length(α)`

# Examples
```jldoctest
julia> B = BracketAlgebra(6,2)
Bracket algebra over P^2 on 6 points with point ordering 1 < 2 < 3 < 4 < 5 < 6 and coefficient ring Integers.

julia> α = [1]
1-element Vector{Int64}:
 1

julia> β = [5,6,2,3]
4-element Vector{Int64}:
 5
 6
 2
 3

julia> γ = [4]
1-element Vector{Int64}:
 4

julia> straightening_sizyge(α, β, γ, B)
[2, 3, 4][1, 5, 6] + [2, 4, 5][1, 3, 6] - [2, 4, 6][1, 3, 5] - [3, 4, 5][1, 2, 6] + [3, 4, 6][1, 2, 5] + [4, 5, 6][1, 2, 3]
```
"""
function straightening_sizyge(α::Vector{<:Integer}, β::Vector{<:Integer}, γ::Vector{<:Integer}, B::BracketAlgebra)
    s = length(α) + 1
    @assert s <= d(B) + 1 "α needs to have length less than or equal to ($(d(B) + 1)), but got $(s-1)."
    @assert length(β) == d(B) + 2 "β needs to have length d + 2 ($(d(B) + 2)), but got $(length(β))."
    @assert length(γ) == d(B) + 1 - s "γ needs to have length d - length(α) ($(d(B) - length(α))), but got $(length(γ))."

    return sum(sign(Perm(vcat(setdiff(collect(1:d(B)+2), τ), τ))) * B(vcat(α, β[setdiff(collect(1:d(B)+2), τ)])) * B(vcat(β[τ], γ)) for τ in combinations(1:d(B)+2, s))
end

"""
    straighten(b::BracketAlgebraElem)

Perform the straightening algorithm on the bracket algebra element `b` as in Stutrmfels 2008, Alg. 3.5.6. The straightening algorithm performs the groebner reduction of the bracket algebra element `b` with the straightening sizyges as a Groebner basis. The result is a normal form of `b` in which every term is standard. See also [`straightening_sizyge`](@ref), [`standard_violation`](@ref), [`is_standard`](@ref).

# Examples
```jldoctest
julia> B = BracketAlgebra(6, 2)
Bracket algebra over P^2 on 6 points with point ordering 1 < 2 < 3 < 4 < 5 < 6 and coefficient ring Integers.

julia> b = B([1,4,5])*B([1,5,6])*B([2,3,4])
[2, 3, 4][1, 5, 6][1, 4, 5]

julia> straighten(b)
[2, 5, 6][1, 4, 5][1, 3, 4] - [3, 5, 6][1, 4, 5][1, 2, 4] + [4, 5, 6][1, 4, 5][1, 2, 3]
```
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

Find all atomic extensors of the bracket algebra elem `b`. These are the equivalence classes of the relation p1 ~ p2 iff substituting p1 = p2 in `b` results in the zero element of the parent bracket algebra.

# Examples
```jldoctest

"""
function atomic_extensors(b::BracketAlgebraElem)
    B = parent(b)

    literals = union(brackets(b)...) # all literals that appear in the bracket expression. Every other literal of B is in its own extensor class, as setting them equal doesn't change the expression.

    extensors = Vector{Int}[[literals[1]]]

    exps = exponent_vectors(b)
    coeffs = coefficients(b)
    function ~(p, q)
        return sum(c * prod(B(replace(B.variables(gens(B.R)[j]), q => p))^e for (j, e) in enumerate(exponent_vector(b, i))) for (i, c) in enumerate(coeffs)) == zero(B)
    end

    # extensors are the equivalence classes of the relation ~. As it is an equivalence relation one only needs to check whether a literal i is equivalent to the first element of every class.
    # The classes also partition the set of literals. Thus, if i is in relation with another literal it doesn't have to be checked against the rest.
    for i in literals[2:end]
        new_class = true # flag that says whether i has to be put into a new equivalence class

        for class in extensors
            if ~(i, class[1])
                push!(class, i)
                new_class = false
                break
            end
        end

        if new_class
            push!(extensors, [i])
        end
    end

    append!(extensors, [[j] for j in setdiff(1:(n(B)), literals)]) # all other literals of B are in their own class.

    return extensors
end