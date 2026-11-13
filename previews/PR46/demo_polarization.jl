using PlasmaBO
using PlasmaBO: c0, ε0, qe, me, mp
using LinearAlgebra
using Printf

println("=================================================================")
println("    PlasmaBO.jl Wave Polarization and Handedness Demo            ")
println("=================================================================")

# 1. Physics Setup
B0 = 1.0e-4          # Background magnetic field along +z [Tesla]
n = 1.0e11           # Plasma density [m^-3]
T = 0.0              # Temperature [eV] - cold limit for fluid

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

println("\n1. Constructing the dispersion matrix using BOFluid solver...")
M = dispersion_matrix(species, B0, kx, kz, BOFluid)

println("2. Solving the eigenvalue problem (eigen!)...")
decomp = eigen!(M)
ωs = decomp.values
V = decomp.vectors

# Analytical Frequencies
wpe2 = n * qe^2 / (ε0 * me)
wpi2 = n * qe^2 / (ε0 * mp)
wce = qe * B0 / me     # Signed (negative for electrons)
wci = qe * B0 / mp     # Signed (positive for protons)
wpe = sqrt(wpe2)
wpi = sqrt(wpi2)

# Dispersion function residuals
D_R(ω) = 1.0 - (c0*k/ω)^2 - wpe2 / (ω * (ω - abs(wce))) - wpi2 / (ω * (ω + wci))
D_L(ω) = 1.0 - (c0*k/ω)^2 - wpe2 / (ω * (ω + abs(wce))) - wpi2 / (ω * (ω - wci))

println("\nIdentified Modes & Physical Verification against Cold-Plasma Theory:")
println("-------------------------------------------------------------------------------------------------------------------------------------------")
@printf("%-5s | %-24s | %-12s | %-10s | %-12s | %-10s | %-10s | %-36s\n", "Index", "Frequency ω (rad/s)", "E-field Norm", "Handedness", "Pol Ratio", "Residual D", "Status", "Physical Description")
println("-------------------------------------------------------------------------------------------------------------------------------------------")

# Sort all 14 solutions by |real(ω)|
sorted_indices = sort(1:length(ωs), by = i -> abs(real(ωs[i])))

for i in sorted_indices
    ω = ωs[i]
    v = V[:, i]
    
    # Extract components and analyze polarization
    Ex, Ey, Ez = electric_field(v)
    Bx, By, Bz = magnetic_field(v)
    E_norm = abs(Ex)^2 + abs(Ey)^2 + abs(Ez)^2
    hand = handedness(v, ω, kx, kz; threshold = 1e-4)
    
    # Residual and status defaults
    residual = NaN
    status = "N/A"
    mode_type = "EM Wave"
    desc = ""
    
    if abs(ω) < 1e1
        mode_type = "Static"
        if abs(Bz) > 0.9
            desc = "Static Longitudinal B-field (Bz)"
        else
            desc = "Static Neutral Plasma Perturbation"
        end
    else
        # Calculate residual to verify if the computed eigenvalue satisfies 
        # the cold-plasma dispersion relation for its physical polarization handedness
        if real(ω) >= 0
            if hand == :R
                residual = abs(D_R(ω))
            elseif hand == :L
                residual = abs(D_L(ω))
            end
        else
            # For negative frequencies, R-mode satisfies D_L(ω) = 0 and L-mode satisfies D_R(ω) = 0
            if hand == :R
                residual = abs(D_L(ω))
            elseif hand == :L
                residual = abs(D_R(ω))
            end
        end
        
        if !isnan(residual)
            status = residual < 1e-4 ? "VERIFIED" : "DISCREPANT"
        end
        
        if E_norm < 1e-16
            mode_type = "Uncoupled Gyro"
            residual = NaN
            status = "N/A"
            if abs(abs(ω) - abs(wci)) < 1e2
                desc = "Proton Cyclotron Gyration (Resonance)"
            else
                desc = "Uncoupled Species Gyration"
            end
        elseif hand == :linear
            mode_type = "Electrostatic"
            desc = "Langmuir Wave " * (real(ω) > 0 ? "(+z)" : "(-z)")
        else
            # EM Waves
            pol = hand == :R ? "R-mode" : "L-mode"
            dir = real(ω) > 0 ? "+z" : "-z"
            if abs(real(ω)) > 2e8
                desc = "EM Light Wave ($(pol), $(dir))"
            elseif abs(real(ω)) > 1e6
                desc = "Whistler Wave ($(pol), $(dir))"
            else
                desc = "Ion Cyclotron Wave ($(pol), $(dir))"
            end
        end
    end
    
    # Format residual string
    res_str = isnan(residual) ? "N/A" : @sprintf("%.2e", residual)
    pol_ratio = polarization_ratio(v, kx, kz)
    # Use absolute value for display; handle non-finite values gracefully
    display_ratio = (isfinite(pol_ratio) ? abs(pol_ratio) : NaN)
    @printf("%-5d | %-24s | %-12.2e | %-10s | %-12.2e | %-10s | %-10s | %-36s\n",
            i,
            @sprintf("%.3e + %.3eim", real(ω), imag(ω)),
            E_norm,
            string(hand),
            display_ratio,
            res_str,
            status,
            desc)
end
println("------------------------------------------------------------------------------------------------------------------")
println("\nSummary of Physics Matching:")
println("  - Electron plasma frequency (ωpe)  : ", @sprintf("%.3e", wpe), " rad/s")
println("  - Ion cyclotron frequency (ωci)    : ", @sprintf("%.3e", wci), " rad/s")
println("  - Electron cyclotron frequency(|ωce|): ", @sprintf("%.3e", abs(wce)), " rad/s")
println("  - High-frequency EM wave speed (c*k) : ", @sprintf("%.3e", c0 * k), " rad/s")
println("\nNote: Handedness is defined relative to the background magnetic field B0 (+z)")
println("      assuming the wave convention e^(i*(k_z * z - ω * t)).")
println("=================================================================")
