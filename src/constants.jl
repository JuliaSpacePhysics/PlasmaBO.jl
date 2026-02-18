# Physical constants
module Constants
using Unitful

export c0, q, mp, me, kb, ε0, μ0, qe

const c0 = ustrip(upreferred(Unitful.c0))

for f in (:μ0, :ε0, :q, :k, :me, :mp)
    @eval const $f = ustrip(upreferred($Unitful.$f))
end

const kb = k
const qe = q
end
