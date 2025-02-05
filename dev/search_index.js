var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = BracketAlgebras","category":"page"},{"location":"#BracketAlgebras","page":"Home","title":"BracketAlgebras","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BracketAlgebras.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BracketAlgebras]","category":"page"},{"location":"#BracketAlgebras.Tabloid-Tuple{BracketAlgebraElem}","page":"Home","title":"BracketAlgebras.Tabloid","text":"tabloid(b::BracketAlgebraElem, ordering::Vector{Int}=collect(1:parent(b).n))\n\nReturn the tabloid corresponding to the bracket monomial b ordered by ordering.  ordering determines the order in which \n\nExample: b = [1,2,3][3,4,5], ordering = [4,3,1,2,5] result: [4 3 5; 3 1 2]\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.atomic_extensors-Tuple{BracketAlgebraElem}","page":"Home","title":"BracketAlgebras.atomic_extensors","text":"atomic_extensors(b::BracketAlgebraElem)\n\nFind all atomic extensors of the bracket algebra elem b. These are the equivalence classes of the relation p1 ~ p2 iff substituting p1 = p2 in b results in the zero element of the bracket algebra.\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.d-Tuple{BracketAlgebra}","page":"Home","title":"BracketAlgebras.d","text":"d(B::BracketAlgebra)\n\nReturn the dimension of the projective space mathbbP^d underlying the bracket algebra B.\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.n-Tuple{BracketAlgebra}","page":"Home","title":"BracketAlgebras.n","text":"n(B::BracketAlgebra)\n\nReturn the number of points considered in the bracket algebra B.\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.point_ordering-Tuple{BracketAlgebra}","page":"Home","title":"BracketAlgebras.point_ordering","text":"point_ordering(B::BracketAlgebra)\n\nReturn the ordering of the n points in the bracket algebra B, that induces the tableaux order.\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.standard_violation-Tuple{Tabloid}","page":"Home","title":"BracketAlgebras.standard_violation","text":"standard_violation(t::Tabloid)\n\nReturn the index of the first violation to the standardness of t, i.e. the first index where t[i,j] > t[i+1, j] with regard to the ordering. Otherwise return nothing. Example: t = [1 2 3; 1 4 5; 1 5 6; 2 3 4] with ordering [1,2,3,4,5,6] => return (3,2) t = [1 2 3; 1 2 4] with ordering [1,2,3,4] => return nothing\n\n\n\n\n\n","category":"method"},{"location":"#BracketAlgebras.straighten-Tuple{BracketAlgebraElem}","page":"Home","title":"BracketAlgebras.straighten","text":"straighten(b::BracketAlgebraElem, ordering::AbstractVector{<:Integer} = collect(1:parent(b).n); chek::Bool = true)\n\nPerform the straightening algorithm to the BracketAlgebraElem b with monomial ordering induced by ordering[1] < ordering[2] < ...\n\nFor details see Sturmfels 2008, chapter 3.1, Handbook of Geometric Constraint System Principles, chapater 4.3\n\n\n\n\n\n","category":"method"}]
}
