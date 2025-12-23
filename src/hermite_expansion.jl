"""
Hermite expansion utilities for arbitrary velocity distributions.

Reference: https://docs.sciml.ai/PolyChaos/stable/orthogonal_polynomials_canonical
"""

"""
    hermite_H(n::Int, x)

Compute physicists' Hermite polynomial H_n(x) using recurrence relation.

# Recurrence Relations
- H_0(x) = 1
- H_1(x) = 2x
- H_{n+1}(x) = 2x*H_n(x) - 2n*H_{n-1}(x)

# Examples
```julia
hermite_H(0, 1.5)  # Returns 1.0
hermite_H(1, 1.5)  # Returns 3.0 (= 2*1.5)
hermite_H(2, 1.5)  # Returns 7.0 (= 2*1.5*3.0 - 2*1*1.0)
```
"""
function hermite_H(n::Int, x::T) where {T <: Number}
    @assert n >= 0 "Hermite polynomial order must be non-negative"

    n == 0 && return one(T)
    n == 1 && return 2 * x

    # Use recurrence relation
    H_prev2 = one(T)
    H_prev1 = 2 * x

    for k in 2:n
        H_curr = 2 * x * H_prev1 - 2 * (k - 1) * H_prev2
        H_prev2 = H_prev1
        H_prev1 = H_curr
    end

    return H_prev1
end


"""
    hermite_coefficients_matrix(nmax)

Compute coefficient matrix for expanding Hermite polynomials in power basis given the maximum order `nmax`.

Returns matrix `cHn[n+1, k+1]` such that:
```
H_n(x) / √(2^n n! √π) = Σ_k cHn[n+1, k+1] * x^k / √(2^{n-k} n! √π)
```

This is used to convert between normalized Hermite basis and power-law basis.

# Returns
- `cHn`: (nmax+1) × (nmax+1) coefficient matrix
"""
function hermite_coefficients_matrix(nmax)
    # Build coefficients for physicists' Hermite polynomials via polynomial recurrence:
    # H_0(x)=1, H_1(x)=2x, H_{n+1}(x)=2x H_n(x) - 2n H_{n-1}(x)
    cHn0 = zeros(Float64, nmax + 1, nmax + 1)
    cHn0[1, 1] = 1.0
    if nmax >= 1
        cHn0[2, 2] = 2.0
    end

    for n in 1:(nmax - 1)
        # Compute H_{n+1} coefficients from H_n and H_{n-1}
        for k in 0:(n + 1)
            term1 = k >= 1 ? 2.0 * cHn0[n + 1, k] : 0.0
            term2 = 2.0 * n * cHn0[n, k + 1]
            cHn0[n + 2, k + 1] = term1 - term2
        end
    end

    # Convert to coefficients consistent with Matlab funa0lm2alm.m
    # Matlab uses a k-dependent normalization (see funa0lm2alm.m lines 39-42):
    #   cHn(n,k) = cHn0(n,k) / sqrt(2^((n-1)-(k-1)) * (n-1)! * sqrt(pi))
    # with 1-based indexing and Hermite order = n-1, power = k-1.
    cHn = zeros(Float64, nmax + 1, nmax + 1)
    for n in 0:nmax
        for k in 0:n
            cHn[n + 1, k + 1] = cHn0[n + 1, k + 1] / sqrt(2.0^(n - k) * factorial(n) * sqrt(π))
        end
    end

    return cHn
end


"""
    hermite_a0_to_a(a0lm)

Convert normalized Hermite basis coefficients: `a0lm` -> alm

Transforms expansion:
```
f(z,x) = Σ_{l,m} a0_{l,m} * ρ_l(z) * u_m(x)
```
to:
```
f(z,x) = Σ_{l,m} a_{l,m} * g_l(z) * h_m(x)
```

where:
- ρ_l, u_m are normalized Hermite functions
- g_l(z) ∝ z^l * exp(-z²/2)
- h_m(x) ∝ x^m * exp(-x²/2)
"""
function hermite_a0_to_a(a0lm)
    lmax, mmax = size(a0lm)
    nmax = max(lmax, mmax)
    cHn = hermite_coefficients_matrix(nmax - 1)
    #   alm = C_z' * a0lm * C_x
    # with C_z = cHn[1:lmax, 1:lmax] and C_x = cHn[1:mmax, 1:mmax].
    Cz = @view cHn[1:lmax, 1:lmax]
    Cx = @view cHn[1:mmax, 1:mmax]
    return @tullio alm[i, j] := Cz[k, i] * a0lm[k, l] * Cx[l, j]
end

@inline function hermite_basis(ξ, l)
    return hermite_H(l, ξ) * exp(-ξ^2 / 2) / sqrt(2^l * factorial(l) * sqrt(π))
end

@inline function hermite_basis(x, x0, μ, l)
    ξ = sqrt(2) * (x - x0) / μ
    return hermite_basis(ξ, l)
end

"""
    hermite_expansion(
        fv::AbstractMatrix{T},
        vz, vx, vtz, vtp;
        Nz = 16, Nx = 16,
        dz = zero(T), dx = zero(T)
    ) where {T}

Compute Hermite expansion coefficients from gridded distribution function.

Expands arbitrary 2D distribution f(v_parallel, v_perp) as:
```
f(vz, vx) = Σ_{l=0}^{Nz} Σ_{m=0}^{Nx} a_{l,m} * ρ_l(vz) * u_m(vx)
```

where ρ_l and u_m are normalized Hermite basis functions. `Nz`, `Nx` are the maximum parallel and perpendicular Hermite indices, respectively.

# Arguments
- `fv`: Distribution function values on (vz, vx) grid
- `vz`: Parallel velocity grid (1D or matches fv size)
- `vx`: Perpendicular velocity grid (1D or matches fv size)
- `vtz`: Parallel thermal velocity (default: 1.0)
- `vtp`: Perpendicular thermal velocity (default: 1.0)
- `dz`: Parallel drift velocity (default: 0.0)
- `dx`: Perpendicular drift velocity (default: 0.0)

# Returns
Named tuple with:
- `alm`: (Nz+1) × (Nx+1) coefficient matrix in power-law basis
- `a0lm`: (Nz+1) × (Nx+1) coefficient matrix in Hermite basis

# Algorithm
1. Compute a0_{l,m} via numerical integration:
   ```
   a0_{l,m} = ∫∫ f(z,x) ρ_l(z) u_m(x) dz dx * (2/(Lz*Lx))
   ```

2. Normalize by distribution integral

3. Convert to power-law basis: a_{l,m} = hermite_a0_to_a(a0_{l,m})
"""
function hermite_expansion(
        fv::AbstractArray{T},
        vz, vx, vtz, vtp; Nz = 16, Nx = 16,
        dz = zero(T), dx = zero(T)
    ) where {T}

    # Accept either 1D vectors (preferred) or prebuilt 2D grids.
    vz_vec = ndims(vz) == 1 ? vz : view(vz, :, 1)
    vx_vec = ndims(vx) == 1 ? vx : view(vx, 1, :)

    @assert size(fv, 1) == length(vz_vec) "Grid size mismatch: fv rows != length(vz)"
    @assert size(fv, 2) == length(vx_vec) "Grid size mismatch: fv cols != length(vx)"

    # Grid spacing (assuming uniform)
    dvz = abs(vz_vec[2] - vz_vec[1])
    dvx = abs(vx_vec[2] - vx_vec[1])

    # Normalized Hermite basis functions
    u_m(x, m) = hermite_basis(x, dx, vtp, m)

    # Precompute basis matrices.
    nx = length(vx_vec)
    Ueff = zeros(T, nx, Nx + 1)
    lranges = 0:Nz
    @tullio ρmat[i, lp1] := hermite_basis(vz_vec[i], dz, vtz, lranges[lp1])

    vx_is_half_space = minimum(vx_vec) >= 0
    vx0_is_first = abs(vx_vec[1]) <= sqrt(eps(real(T))) * max(one(T), abs(vx_vec[end]))

    @views if vx_is_half_space
        # Matlab: treat the first column (usually vx=0) separately.
        if vx0_is_first
            @tullio Ueff[i, mp] = (
                i == 1 ? u_m(vx_vec[1], mp - 1) :
                    (u_m(vx_vec[i], mp - 1) + u_m(-vx_vec[i], mp - 1))
            )
        else
            @tullio Ueff[i, mp] = u_m(vx_vec[i], mp - 1) + u_m(-vx_vec[i], mp - 1)
        end
    else
        @tullio Ueff[i, mp] = u_m(vx_vec[i], mp - 1)
    end

    a0lm = ρmat' * fv * Ueff
    a0lm .*= (dvz * dvx * 2 / (vtz * vtp))

    # Normalization factor
    bs = dx / vtp
    As = exp(-bs^2) + sqrt(π) * bs * erfc(-bs)
    cs0 = 1 / (sqrt(π^3) * vtz * vtp^2 * As)

    # Convert to power-law basis
    alm = hermite_a0_to_a(a0lm)
    alm ./= cs0

    return (;
        alm = alm,
        a0lm = a0lm,
    )
end
