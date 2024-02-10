# Test if derivatives are correct
using FrankWolfe
using FiniteDifferences
using LinearAlgebra
using Random
using ODWB
using Test

"""
Check if the gradient using finite differences matches the grad! provided.
Copied from FrankWolfe package: https://github.com/ZIB-IOL/FrankWolfe.jl/blob/master/examples/plot_utils.jl
"""
function check_gradients(grad!, f, gradient, num_tests=10, tolerance=1.0e-5)
    for i in 1:num_tests
        random_point = rand(length(gradient))
        grad!(gradient, random_point)
        if norm(grad(central_fdm(5, 1), f, random_point)[1] - gradient) > tolerance
            @warn "There is a noticeable difference between the gradient provided and
            the gradient computed using finite differences.:\n$(norm(grad(central_fdm(5, 1), f, random_point)[1] - gradient))"
            return false
        end
    end
    return true
end


@testset "Derivative A-opt" begin
    for dim in [20,50,80]
        n = Int(floor(dim/4))
        @show dim, n
        gradient = rand(dim)
        A, _, _, _ = ODWB.build_data(seed,dim, n, false, false)
        f, grad! = ODWB.build_a_criterion(A, false, μ=1e-2)

        @test check_gradients(grad!, f, gradient)
    end
end

@testset "Derivative D-opt" begin
    for dim in [20,50,80]
        n = Int(floor(dim/4))
        @show dim, n
        gradient = rand(dim)
        A, _, _, _ = ODWB.build_data(seed,dim, n, false, false)
        f, grad! = ODWB.build_d_criterion(A, false)

        @test check_gradients(grad!, f, gradient)
    end
end

@testset "Derivative A-Fusion" begin
    for dim in [20,50,80]
        n = Int(floor(dim/4))
        @show dim, n
        gradient = rand(dim)
        A, C, _, _ = ODWB.build_data(seed,dim, n, true, false)
        f, grad! = ODWB.build_a_criterion(A, true, C=C)

        @test check_gradients(grad!, f, gradient)
    end
end

@testset "Derivative D-Fusion" begin
    for dim in [20,50,80]
        n = Int(floor(dim/4))
        @show dim, n
        gradient = rand(dim)
        A, C, _, _ = ODWB.build_data(seed,dim, n, true, false)
        f, grad! = ODWB.build_d_criterion(A,true, C=C)

        @test check_gradients(grad!, f, gradient)
    end
end

@testset "Derivatives Custom BB solver" begin
    for p in 0:-0.5:-2
        @show p
        for dim in [20,50,80]
            n = Int(floor(dim/4))
            @show dim, n
            gradient = rand(dim)
            A, _, _, _ = ODWB.build_data(seed,dim, n, true, false)
            f, grad! = ODWB.build_matrix_means_objective(A,p)

            @test check_gradients(grad!, f, gradient)
        end
    end
end
