Vaje01

Avtor
[Urh Perčič]

Opis
Projekt Vaje01 vsebuje implementacijo algoritma Successive Over-Relaxation (SOR) za razpršene matrike. Cilj je reševanje tridiagonalnega sistema Ax=b s pomočjo SOR iteracije. Projekt vključuje tudi poseben podatkovni tip RazprsenaMatrika, namenjen učinkovitemu shranjevanju in obdelavi razpršenih matrik.

Struktura
src/Vaje01.jl: Glavna implementacija, vključno s podatkovnimi tipi in funkcijami.
test/runtests.jl: Testne datoteke za preverjanje pravilnosti algoritmov.

Uporaba kode
Za uporabo kode, najprej vključite modul s pomočjo:
using .Vaje01

Nato lahko uporabite funkcije, kot je sor za reševanje sistema Ax=b:
A_sparse = SparseMatrixCSC([4.0 -1.0; -1.0 4.0])
A_scattered = Vaje01.ScatteredArray(A_sparse)
b = [0.0, 0.0]
x0 = [0.0, 0.0]
omega = 1.0
tol = 1e-10

optimal_omega, omegas, iterations = Vaje01.find_optimal_omega(A_scattered, b, x0, tol)

Poganjanje testov
Za zagon testov odprite Julia in izvedite
include("test/runtests.jl")
