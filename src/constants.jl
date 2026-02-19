# Physical constants
module Constants
using Unitful

export c0, q, mp, me, kb, ε0, μ0, qe

const c0 = Unitful.c0.val
const ε0 = uconvert(Unitful.F / Unitful.m, Unitful.ε0).val

for f in (:μ0, :q, :k, :me, :mp)
    @eval const $f = ($Unitful.$f).val
end

const kb = k
const qe = q
end
