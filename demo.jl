include("./src/Vaje01.jl")
using .Vaje01
using SparseArrays

function main()
    A_sparse = SparseMatrixCSC([4.0 -1.0; -1.0 4.0])
    A_scattered = Vaje01.ScatteredArray(A_sparse)

    omega = 1.0
    b = [1.0, 1.0]
    x0 = [0.0, 0.0]
    tol = 1e-10

    optimal_omega, omegas, iterations = Vaje01.find_optimal_omega(A_scattered, b, x0, tol)
    omegas_vector = collect(omegas)  
    iterations_vector = collect(iterations) 

    Vaje01.plot_omega_dependency(omegas_vector, iterations_vector)

    println("Finished executing main function.")
end

main()