
# General CRNT functions

export product_matrix, 
reaction_pairs, 
stoichiometric_and_kinetic_matrix_from_reaction_pairs,
minimal_siphons,
critical_minimal_siphons,
has_critical_siphon,
reconstruct_stoichiometric_matrix_with_row_space_and_conserved_quantities


function product_matrix(N::QQMatrix, M::ZZMatrix)
    @req size(N) == size(M) "Stoichiometric and kinetic matrix must have the same size"
    @req all(is_one, denominator.(N)) "Stoichiometric matrix needs to have integer entries"
    N = ZZ.(N)
    return M + N
end

function reaction_pairs(N::QQMatrix, M::ZZMatrix)
    P = product_matrix(N, M)
    reactions = [(reactant=M[:, j], product=P[:, j]) for j in 1:size(M, 2)]
    return reactions
end

function stoichiometric_and_kinetic_matrix_from_reaction_pairs(reactions::Vector{NamedTuple{(:reactant, :product),Tuple{Vector{T},Vector{T}}}}) where {T}
    @assert all(length(reactions[1].reactant) == length(rxn.reactant) && length(reactions[1].reactant) == length(rxn.product) for rxn in reactions) "All vectors in all pairs must have the same length"
    @assert T == ZZRingElem || T == Int "Entries must be of type ZZRingElem or Int"
    M = matrix(ZZ, hcat([rxn.reactant for rxn in reactions]...))
    N = matrix(QQ, hcat([rxn.product - rxn.reactant for rxn in reactions]...))
    return N, M
end



"""
    minimal_siphons(N::QQMatrix, M::ZZMatrix)

Computes the inclusion-minimal siphons of the network with stoichiometric matrix `N` and kinetic matrix `M`. 

"""
function minimal_siphons(N::QQMatrix, M::ZZMatrix)
    n = nrows(M)
    reactions = reaction_pairs(N, M)

    # Helper function to check if Z is a siphon
    function is_siphon(Z)
        if isempty(Z)
            return false  # Exclude the empty set
        end
        for rxn in reactions
            for i in Z
                if rxn.product[i] > 0
                    if !any(j -> rxn.reactant[j] > 0, Z)
                        return false
                    end
                end
            end
        end
        return true
    end

    # Helper function for growing siphons
    function extend_to_siphon(candidate, remaining_indices, list_of_siphons)
        # Add current candidate to the list if it is a siphon and does not contain a smaller siphon
        if is_siphon(candidate)
            if !any(s -> issubset(s, candidate) && s != candidate, list_of_siphons)
                push!(list_of_siphons, candidate)
            end
        end
        for i in remaining_indices
            new_candidate = union(candidate, [i])
            new_remaining = filter(x -> x > i, remaining_indices)
            extend_to_siphon(new_candidate, new_remaining, list_of_siphons)
        end
        return list_of_siphons 
    end

    return extend_to_siphon(Set{Int}(), 1:n, Set{Int}[])
end


"""
    has_critical_siphon(N::QQMatrix, M::ZZMatrix)

Checks whether a network has a critical siphon (i.e., a siphon that does not contain 
the support of a positive conservation law).

If there is a critical minimal siphon, the function returns `true`.

If there are no critical minimal siphons, the function returns `false`. 
In this case, the network lacks relevant boundary steady states and displays persistence.

"""
function has_critical_siphon(N::QQMatrix, M::ZZMatrix)
    return !isempty(critical_minimal_siphons(N, M))
end


"""
    critical_minimal_siphons(N::QQMatrix, M::ZZMatrix)

Find the critical minimal siphons of the network with stoichiometric matrix `N` and kinetic matrix `M`
(i.e., the minimal siphons that do not contain the support of a positive conservation law).

"""
function critical_minimal_siphons(N::QQMatrix, M::ZZMatrix)
    siphons = minimal_siphons(N, M)
    W = kernel(N, side=:left)
    critical_siphons = Vector{Set{Int}}()
    for Z in siphons
        # Check if there is a nonnegative vector in rowspan(W) with support contained in Z
        if !_has_nonnegative_rowspace_vector_supported_in(W, collect(Z))
            push!(critical_siphons, Z)
        end
    end
    return critical_siphons
end

# Helper function: Checks whether A has a nonzero nonnegative vector x in its row space whose support is contained in Z
function _has_nonnegative_rowspace_vector_supported_in(W::Union{ZZMatrix, QQMatrix}, Z::Vector{Int})
    @req !isempty(Z) "Z must be nonempty"
    @req all(i -> 1 <= i <= ncols(W), Z) "Z must be a subset of 1:ncols(W)"

    C = kernel(W, side=:right)

    # Inequalities: impose x_Z >= 0
    inequalities = (-identity_matrix(QQ, length(Z)), zeros(Int, length(Z)))

    # Normalization: impose sum(x_Z) = 1 (ruling out the zero vector)
    normalization = matrix(QQ, 1, length(Z), ones(Int, length(Z)))

    # Equalities: impose that x obtained by extending x_Z by zeros lies in rowspan(A)
    equalities = (
        vcat(transpose(C[Z, :]), normalization),
        vcat(zeros(Int, ncols(C)), [1])
    )

    # Check for feasibility
    P = polyhedron(inequalities, equalities)
    return is_feasible(P)
end


function reconstruct_stoichiometric_matrix_with_row_space_and_conserved_quantities(C::QQMatrix, W::QQMatrix)
    @req nrows(W) == rank(W) "Matrix of conserved quantities needs to have full row rank"
    @req nrows(C) == rank(C) "Coefficient matrix needs to have full row rank"
    @req nrows(C) + nrows(W) == ncols(W) "Not enough conserved quantities"
    W = rref(W)[2]
    C = rref(C)[2]
    pivot_colums = [findfirst(row .!= 0) for row in eachrow(W)]
    non_pivot_columns = setdiff(1:ncols(W), pivot_colums)
    N = zero_matrix(QQ,ncols(W),ncols(C))
    N[non_pivot_columns,:] = C 
    N[pivot_colums,:] = -W[:,non_pivot_columns]*C
    # Todo: Remove these checks
    @req row_space(kernel(N, side=:left)) == row_space(W) "Reconstruction failed"
    @req row_space(N) == row_space(C) "Reconstruction failed"
    return N
end
