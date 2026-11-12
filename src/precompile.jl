using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    p, e = Proton(), Electron()
    B0, n, T = 5.0e-9, 5.0e6, 12.94
    maxwellian = Plasma(Population(p, BiMaxwellian; T, n), Population(e, BiMaxwellian; T, n); B0)
    at = Wavevector(; k = 1.0e-4, theta = 0.5)
    sweep = AngleSweep(; k = [1.0e-4, 2.0e-4], theta = 0.5)
    @compile_workload begin
        solve(DispersionProblem(maxwellian, sweep), BOHH())
        solve(DispersionProblem(maxwellian, Seed(0.1 + 0.01im, at), sweep), BOHH())
        dispersion_matrix(maxwellian, at, BOHH())
        HermiteVDF(ProductBiKappa(p; T_perp = T, T_para = T, kappa_perp = 5.5, kappa_para = 5.5))
    end
end
