include("../src/Vaje01.jl")
using Test
using .Vaje01
using SparseArrays
using LinearAlgebra 
using Coverage

@testset "ScatteredArray Tests" begin
    @testset "Construction and Indexing" begin
        A_sparse = SparseMatrixCSC([1.0 2.0; 3.0 4.0])
        A_scattered = ScatteredArray(A_sparse)

        @test A_scattered[1, 1] == 1.0
        @test A_scattered[1, 2] == 2.0
        @test A_scattered[2, 1] == 3.0
        @test A_scattered[2, 2] == 4.0
        @test A_scattered[1, 3] == 0.0  
    end

    @testset "SOR Method Tests" begin
        A_sparse = SparseMatrixCSC([4.0 -1.0; -1.0 4.0])
        A_scattered = ScatteredArray(A_sparse)
        b = [5.0, 11.0]
        x0 = [0.0, 0.0]
        omega = 1.0

        result, _ = sor(A_scattered, b, x0, omega)
        known_solution = A_sparse \ b
        @test norm(result - known_solution) < 1e-4
    end

    @testset "Omega Value Variation" begin
        omegas = [0.5, 1.0, 1.5]
        A_sparse = SparseMatrixCSC([1.0 2.0; 3.0 4.0])
        A_scattered = ScatteredArray(A_sparse)
        b = [5.0, 11.0]
        x0 = [0.0, 0.0]
    
        for omega in omegas
            @test all(diag(A_sparse) .!= 0)  
            result, _ = sor(A_scattered, b, x0, omega)
            print(result)  
            @test all(isfinite, result)  
        end
    end
end

coverage_results = Coverage.process_folder("../src/Vaje01.jl")
for file in keys(coverage_results)
    println("Coverage for $file")
    for (num, cnt) in coverage_results[file]
        println("Line $num: executed $cnt times")
    end
end
