export plasma_frequency, Debye_length

plasma_frequency(q, n, mass) = abs(q) * sqrt(n / mass / ε0)

gyrofrequency(q, B, mass) = q * B / mass

const wp_ = plasma_frequency


Debye_length(n, T) = sqrt(ε0 * kb * T / (n * q^2))

Debye_length(s::Species) = Debye_length(s.n, s.Tz * q / kb)

Debye_length(species) = sqrt(1.0 / sum(1.0 ./ Debye_length.(species) .^ 2))
