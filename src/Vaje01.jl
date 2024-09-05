module Vaje01
    using LinearAlgebra
    using SparseArrays
    using Plots

    struct ScatteredArray{T <: Number}
        V::Vector{Vector{T}}
        I::Vector{Vector{Int}}
    end

    function ScatteredArray(A::SparseMatrixCSC{T}) where T <: Number
        V = [Vector{T}() for _ in 1:size(A, 1)]
        I = [Vector{Int}() for _ in 1:size(A, 1)]
        for col in 1:size(A, 2)
            for row_idx in nzrange(A, col)
                row = rowvals(A)[row_idx]
                val = nonzeros(A)[row_idx]
                push!(V[row], val)
                push!(I[row], col)
            end
        end
        return ScatteredArray{T}(V, I)
    end

    function Base.getindex(A::ScatteredArray{T}, i::Int, j::Int) where T
        idx = findfirst(==(j), A.I[i])
        return idx !== nothing ? A.V[i][idx] : zero(T)
    end

    function Base.setindex!(A::ScatteredArray{T}, value, i::Int, j::Int) where T
        idx = findfirst(==(j), A.I[i])
        if value != 0
            if idx !== nothing
                A.V[i][idx] = value
            else
                push!(A.V[i], value)
                push!(A.I[i], j)
            end
        else
            if idx !== nothing
                deleteat!(A.V[i], idx)
                deleteat!(A.I[i], idx)
            end
        end
    end

    Base.firstindex(A::ScatteredArray) = 1
    Base.lastindex(A::ScatteredArray{T}) where T = length(A.V)

    function Base.:*(A::ScatteredArray{T}, x::Vector{T}) where T
        n = length(A.V)
        result = zeros(T, n)
        for i in 1:n
            for k in 1:length(A.I[i])
                result[i] += A.V[i][k] * x[A.I[i][k]]
            end
        end
        return result
    end

    function sor(A::ScatteredArray{T}, b::Vector{T}, x0::Vector{T}, omega::Float64, tol::Float64=1e-10) where T
        x = copy(x0)
        n = length(b)
        for iter in 1:10000
            x_old = copy(x)
            for i in 1:n
                sigma = sum(A[i, j] * x[j] for j in 1:n if j != i)
                if A[i, i] == 0
                    error("Zero diagonal element at index $i, which will cause division by zero.")
                end
                x[i] = (1 - omega) * x_old[i] + omega * (b[i] - sigma) / A[i, i]
            end
            if norm(x - x_old, Inf) < tol
                return x, iter
            end
        end
        return x, 10000
    end
    
    function find_optimal_omega(A::ScatteredArray, b::Vector, x0::Vector, tol::Float64)
        println("Searching for optimal omega...")
        optimal_omega = 0.0
        min_iterations = Inf
        omega_values = 0.1:0.1:1.9
        iterations = Int64[]  
    
        for omega in omega_values
            _, it = sor(A, b, x0, omega, tol)
            push!(iterations, it) 
            if it < min_iterations
                min_iterations = it
                optimal_omega = omega
            end
        end
    
        return optimal_omega, omega_values, iterations
    end

    function plot_omega_dependency(omegas::AbstractVector, iterations::Vector{Int})
        plot_object = plot(omegas, iterations, xlabel="Omega", ylabel="Iterations", title="SOR Convergence")
        display(plot_object)
    end

    export ScatteredArray, sor, find_optimal_omega, plot_omega_dependency
end


