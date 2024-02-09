# Test for strong convexity
using Random
using LinearAlgebra
using LinearMaps
using Kronecker

# https://github.com/MichielStock/Kronecker.jl/issues/91
#khatri_rao(A::AbstractMatrix, B::AbstractMatrix) = hcat(map(LinearMap∘kronecker, eachcol(A), eachcol(B))...)
khatri_rao(A, B) = mapreduce(kron, hcat, eachcol(A), eachcol(B))

#khatri_rao(A::AbstractMatrix, B::AbstractMatrix) = hcat(map(⊗, eachcol(A), eachcol(B))...)

function column_wise_kronecker(A,B)
    m_A, n_A = size(A)
    m_B, n_B = size(B)

    @assert n_A == n_B

    C = zeros(m_B*m_A, n_A)
    for i in 1:n_A
        H = A[:,i] * B[:,i]'
        for j in 1:m_A
            C[m_B*(j-1)+1:j*m_B,i] = H[j,:]
        end
    end
    return C
end

dimensions = [20,30,40,50,60,80,100,120,150,200] #[30,40]
n_set = [4,5,6,8,10]
maxseed = 5
full_rank_count = 0

# Rank if adding a constant to all entries
for m in dimensions
    for k in n_set
        local n = Int(round(m/k))
        @show m, n
        for seed in 1:10
            Random.seed!(seed)

            M = randn(m,n)

            @assert rank(M) == n

            local c = abs(minimum(M)) + 1

            N = M .+ c

            if rank(N) == n
               global full_rank_count += 1
            end
            @show seed, rank(M), rank(N)
        end 
    end
end

@show 10*length(dimensions)*length(n_set), full_rank_count, full_rank_count/(10*length(dimensions)*length(n_set))

println("\n")

full_rank_count = 0

# Rank with Kronecker product
for m in dimensions
    for k in n_set
        count = 0
        local n = Int(round(m/k))
        println("Dimesions m=$(m) and n=$(n)")
        println("n² = $(n^2) ≥ m=$(m): $(n^2 ≥ m)")
        println("n choose 2 = $(n^2-n) ≥ m =$(m): $(n^2-n ≥ m)")
        for seed in 1:maxseed
            Random.seed!(seed)
            local A = rand(m,n)
            C = column_wise_kronecker(A',A')
            if rank(C) == m
                global full_rank_count += 1
                count += 1
            end
            @show seed, rank(A), rank(C)
        end
        @show maxseed, count, count/maxseed
    end
    println("\n")
end

println("Total number: $(length(dimensions)*length(n_set)*maxseed)  Full Rank Count: $(full_rank_count) Ratio: $(full_rank_count/(length(dimensions)*length(n_set)*maxseed))")