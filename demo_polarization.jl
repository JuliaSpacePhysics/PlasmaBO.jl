using PlasmaBO
using LinearAlgebra

println("=================================================================")
println("    PlasmaBO.jl Wave Polarization and Handedness Demo            ")
println("=================================================================")

eigenvectors = true

# 1. Physics Setup
B0 = 1.0e-4          # Background magnetic field along +z [Tesla]
n = 1.0e11           # Plasma density [m^-3]
T = 1.0              # Temperature [eV]

# Parallel propagation (along +z direction)
k = 1.0
kx = 0.0
kz = k

println("Parameters:")
println("  - Background Magnetic Field (B0): ", B0, " T (along +z)")
println("  - Plasma Density (n)            : ", n, " m^-3")
println("  - Temperature (T)               : ", T, " eV")
println("  - Wave Vector (kx, kz)          : (", kx, ", ", kz, ") m^-1")

# Define species
e_vdf = Maxwellian(:e, n, T)
species = [e_vdf, FluidSpecies(:p, n, T)]

println("\nSolving for waves using BOFluid (electromagnetic fluid solver)...")
decomp = solve(species, B0, kx, kz, BOFluid, eigenvectors)

ωs = decomp.values
V = decomp.vectors

println("\nIdentified Modes & Their Polarizations:")
println("---------------------------------------------------------------------------------")
printf_fmt = "%-5s | %-32s | %-12s | %-12s | %-10s\n"
@eval import Printf: @printf
@printf("%-5s | %-32s | %-24s | %-12s\n", "Index", "Frequency ω (rad/s)", "Polarization Ratio (Ey/Ex)", "Handedness")
println("---------------------------------------------------------------------------------")

for i in 1:length(ωs)
    ω = ωs[i]
    v = V[:, i]
    
    # Extract components
    Ex, Ey, Ez = electric_field(v)
    P = polarization_ratio(v)
    hand = handedness(v; threshold = 1e-4)
    
    # We filter out zero/near-zero frequencies for readability
    if abs(ω) > 1e1
        @printf("%-5d | %-32s | %-24s | %-12s\n", 
                i, 
                string(round(real(ω), digits=3)) * " + " * string(round(imag(ω), digits=3)) * "im",
                string(round(real(P), digits=3)) * " + " * string(round(imag(P), digits=3)) * "im",
                string(hand))
    end
end
println("---------------------------------------------------------------------------------")
println("\nNote: Handedness is defined relative to the background magnetic field B0 (+z).")
println("=================================================================")
