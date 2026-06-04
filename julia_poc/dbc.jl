# ============================================================================
# dbc.jl  —  Pure-Julia SZ-DBC detailed-balance checker (Phase 1 BFS + DB check)
# ============================================================================
# This is the Julia-native implementation of the SZ-DBC pipeline. When invoked
# via check.jl on a translated algorithm (j_examples/*.jl) it runs the ENTIRE
# checker in Julia with no Mathematica dependency:
#
#     state enumeration -> tau-augmented BFS (translational-invariance check +
#     exact symbolic path enumeration) -> ergodicity reachability ->
#     detailed-balance check (HiGHS LP chamber enumeration + exact rational
#     grouping).
#
# DESIGN PRINCIPLE: correctness over speed. The checker must NEVER return a
# false PASS (DB satisfied when it is not) or a false FAIL. Anything the Julia
# engine cannot represent EXACTLY is raised as a `CantHandle` error and aborts
# the run loudly, rather than producing a possibly-wrong verdict.
#
# See README.md ("Julia-native checker") for the supported-primitive list and
# the documented edge cases that trigger a hard error.
# ============================================================================

using Printf

# ----------------------------------------------------------------------------
# Exceptions
# ----------------------------------------------------------------------------
struct OutOfBitsException  <: Exception end          # BFS path needs more bits
struct OutOfRangeException <: Exception end          # rejected bit pattern (rejection sampling)
struct CantHandle          <: Exception; msg::String end   # unsupported -> hard abort

cant(msg) = throw(CantHandle(msg))

# ============================================================================
# SECTION 1 — TauNum: exact linear tau-tracking with nonlinear taint
# ============================================================================
# Represents a scalar of the form
#
#     v + cr*tau_r + cc*tau_c        (+ tainted higher-order tau terms)
#
# where v, cr, cc are exact rationals and `tainted` records that a genuinely
# nonlinear tau term (tau_r^2, tau_r*tau_c, ...) was produced and dropped from
# the linear representation. A value is translation-invariant ("tau-free") iff
#
#     cr == 0 && cc == 0 && !tainted
#
# The tau symbols are a pure detection device: the real algorithm runs at
# tau = 0. After detection, `tau0(x)` extracts the genuine value v (the value
# the real lattice algorithm would compute).
#
# Subtraction of two positions carrying the SAME offset cancels tau exactly
# (cr_a - cr_b = 0), which is how pairwise differences become tau-free.
# Squaring a tau-augmented coordinate (cr != 0) produces a tau_r^2 term and
# sets `tainted`, which is how absolute-position energies (quadratic_field) are
# detected as translation-NON-invariant.
# ----------------------------------------------------------------------------

const Q = Rational{BigInt}

struct TauNum
    v::Q
    cr::Q
    cc::Q
    tainted::Bool
end

TauNum(v::Integer)            = TauNum(Q(v), Q(0), Q(0), false)
TauNum(v::Rational)           = TauNum(Q(v), Q(0), Q(0), false)
tau_const(v)                  = TauNum(Q(v), Q(0), Q(0), false)
# tau-augment a coordinate: value r, unit offset along the chosen axis.
tau_r_aug(r::Integer) = TauNum(Q(r), Q(1), Q(0), false)
tau_c_aug(r::Integer) = TauNum(Q(r), Q(0), Q(1), false)

is_tau_free(x::TauNum) = x.cr == 0 && x.cc == 0 && !x.tainted
tau0(x::TauNum)        = x.v          # value at tau = 0 (the real value)

Base.convert(::Type{TauNum}, v::Integer)  = TauNum(v)
Base.convert(::Type{TauNum}, v::Rational) = TauNum(v)
Base.promote_rule(::Type{TauNum}, ::Type{<:Integer})  = TauNum
Base.promote_rule(::Type{TauNum}, ::Type{<:Rational}) = TauNum

Base.:+(a::TauNum, b::TauNum) = TauNum(a.v+b.v, a.cr+b.cr, a.cc+b.cc, a.tainted|b.tainted)
Base.:-(a::TauNum, b::TauNum) = TauNum(a.v-b.v, a.cr-b.cr, a.cc-b.cc, a.tainted|b.tainted)
Base.:-(a::TauNum)            = TauNum(-a.v, -a.cr, -a.cc, a.tainted)

function Base.:*(a::TauNum, b::TauNum)
    # Linear part of the product; flag any nonlinear tau cross-term as taint.
    v  = a.v*b.v
    cr = a.v*b.cr + a.cr*b.v
    cc = a.v*b.cc + a.cc*b.v
    nonlinear = a.tainted | b.tainted |
                (a.cr != 0 && b.cr != 0) |   # tau_r^2
                (a.cc != 0 && b.cc != 0) |   # tau_c^2
                (a.cr != 0 && b.cc != 0) |   # tau_r*tau_c
                (a.cc != 0 && b.cr != 0)
    TauNum(v, cr, cc, nonlinear)
end

Base.:+(a::TauNum, b::Union{Integer,Rational}) = a + TauNum(b)
Base.:+(a::Union{Integer,Rational}, b::TauNum) = TauNum(a) + b
Base.:-(a::TauNum, b::Union{Integer,Rational}) = a - TauNum(b)
Base.:-(a::Union{Integer,Rational}, b::TauNum) = TauNum(a) - b
Base.:*(a::TauNum, b::Union{Integer,Rational}) = a * TauNum(b)
Base.:*(a::Union{Integer,Rational}, b::TauNum) = TauNum(a) * b

function Base.:^(a::TauNum, n::Integer)
    n < 0 && cant("TauNum raised to a negative power")
    r = TauNum(1)
    for _ in 1:n; r = r * a; end
    r
end

# Equality is exact and tau-aware (used only for canonical dedup of tau-free
# values; comparisons that BRANCH go through require_tau_free first).
Base.:(==)(a::TauNum, b::TauNum) = a.v==b.v && a.cr==b.cr && a.cc==b.cc && a.tainted==b.tainted
Base.hash(a::TauNum, h::UInt) = hash((a.v,a.cr,a.cc,a.tainted), h)

# ============================================================================
# SECTION 2 — Atoms, linear forms, particles, geometry
# ============================================================================

# An Atom is a symbolic parameter the energy / weights are linear in:
#   * a coupling atom  couplingJ[a,b,d2]  (canonicalised so a <= b), or
#   * an "extra" field-like parameter (e.g. fieldH), keyed by a Symbol.
struct Atom
    iscoupling::Bool
    a::Int
    b::Int
    d2::Int
    name::Symbol
end
Jc(a::Int, b::Int, d2::Int) = a <= b ? Atom(true, a, b, d2, :_) : Atom(true, b, a, d2, :_)
Xparam(name::Symbol)        = Atom(false, 0, 0, 0, name)

# Deterministic total order so the global atom list has a stable index.
function Base.isless(x::Atom, y::Atom)
    x.iscoupling != y.iscoupling && return x.iscoupling   # couplings before extras
    if x.iscoupling
        return (x.a, x.b, x.d2) < (y.a, y.b, y.d2)
    else
        return x.name < y.name
    end
end

# A LinForm maps atoms to (possibly tau-dependent) coefficients.
const LinForm  = Dict{Atom, TauNum}      # used while building energies (carries tau)
const RatForm  = Dict{Atom, Q}           # tau=0 form stored in conditions / exponents

function addcoef!(lf::LinForm, a::Atom, c::TauNum)
    nc = get(lf, a, TauNum(0)) + c
    if nc.v == 0 && nc.cr == 0 && nc.cc == 0 && !nc.tainted
        delete!(lf, a)
    else
        lf[a] = nc
    end
    lf
end

function linsub(p::LinForm, q::LinForm)::LinForm     # p - q
    r = copy(p)
    for (a, c) in q; addcoef!(r, a, -c); end
    r
end

any_tau_dep(lf::LinForm) = any(!is_tau_free(c) for c in values(lf))

# Substitute tau = 0 and drop coefficients that vanish at tau = 0.
function tau0_form(lf::LinForm)::RatForm
    r = RatForm()
    for (a, c) in lf
        c.v != 0 && (r[a] = c.v)
    end
    r
end

ratscale(f::RatForm, s::Q)::RatForm = RatForm(a => v*s for (a, v) in f)

# ---- particles / states ----
struct Particle
    r::TauNum
    c::TauNum
    t::Int
end
const PState = Vector{Particle}

# tau-augment a concrete (r,c,type) seed particle along both axes.
aug_particle(r::Int, c::Int, t::Int) = Particle(tau_r_aug(r), tau_c_aug(c), t)

# ---- geometry (all flag tau-violation if used on a tau-dependent value) ----

# Periodic Mod of a coordinate value into 0..n-1, flagging tau dependence.
function pmod(x::TauNum, n::Int)::Int
    is_tau_free(x) || tau_violation!("Mod / branch on an absolute (tau-dependent) position")
    v = tau0(x)
    denominator(v) == 1 || cant("non-integer position coordinate: $v")
    mod(Int(numerator(v)), n)
end

# Squared minimum-image distance between two tau-augmented positions.
# Subtraction cancels tau for genuine pairwise differences; if it does not
# (an absolute-position energy), pmod flags the tau-violation.
function pbc_d2(p::Particle, q::Particle, n::Int)::Int
    dra = pmod(p.r - q.r, n)
    dca = pmod(p.c - q.c, n)
    min(dra, n - dra)^2 + min(dca, n - dca)^2
end

# Are two positions the same site under PBC?  (occupancy test; tau-free diff)
same_site(p::Particle, q::Particle, n::Int) =
    pmod(p.r - q.r, n) == 0 && pmod(p.c - q.c, n) == 0

# ============================================================================
# SECTION 3 — BitSeqRNG and the symbolic random primitives
# ============================================================================
# Mirrors RunWithBitsAT in dbc_core.wl. Each BFS path replays the algorithm
# against a fixed bit string. Selection primitives consume bits and multiply
# the exact rational path coefficient; the Metropolis acceptance records a
# SYMBOLIC factor (the clamped Boltzmann threshold) rather than a float, so the
# downstream DB check is algebraically exact.

# Tau-violation is reported through a process-global flag, reset per BFS run,
# so that the geometry helpers above (which have no rng handle) can raise it.
const _TAU     = Ref(false)
const _TAU_MSG = Ref("")
tau_violation!(msg::String) = (_TAU[] = true; _TAU_MSG[] = msg; nothing)

# One Metropolis accept/reject factor:
#   threshold = Piecewise[{{1, cond_lhs <= 0}}, exp(-beta * expo)]
#   accept branch contributes the threshold; reject branch 1 - threshold.
struct MFactor
    cond::RatForm     # cond_lhs (the clamp condition is cond_lhs <= 0)
    expo::RatForm     # Boltzmann exponent g
    accepted::Bool
end

mutable struct BitSeqRNG
    bits::Vector{Int}
    pos::Int
    coeff::Q                  # rational selection-probability product
    factors::Vector{MFactor}  # symbolic acceptance factors, in order
end
BitSeqRNG(bits::Vector{Int}) = BitSeqRNG(bits, 0, Q(1), MFactor[])

function read_bit!(rng::BitSeqRNG)::Int
    rng.pos += 1
    rng.pos > length(rng.bits) && throw(OutOfBitsException())
    rng.coeff *= 1 // 2
    rng.bits[rng.pos]
end

function read_bits_int!(rng::BitSeqRNG, k::Int)::Int
    acc = 0
    for _ in 1:k
        acc = acc * 2 + read_bit!(rng)
    end
    acc
end

# Bits needed to index 0..n-1, matching Mathematica IntegerLength[n-1, 2].
nbits(n::Int) = n <= 1 ? 0 : ndigits(n - 1; base = 2)

# Uniform choice of an index 1..n with exact rejection-sampling weight 1/n.
function rand_choice_index!(rng::BitSeqRNG, n::Int)::Int
    n == 0 && cant("RandomChoice over an empty list")
    n == 1 && return 1
    k   = nbits(n)
    val = read_bits_int!(rng, k)
    val >= n && throw(OutOfRangeException())
    rng.coeff *= Q(2)^k // n          # undo the (1/2)^k from read_bit!, leave 1/n
    val + 1
end

rand_choice!(rng::BitSeqRNG, list::AbstractVector) = list[rand_choice_index!(rng, length(list))]

# RandomInteger[{lo, hi}] with exact rejection sampling.
function rand_integer!(rng::BitSeqRNG, lo::Int, hi::Int)::Int
    n = hi - lo + 1
    n <= 0 && cant("RandomInteger[{$lo,$hi}]: inverted range (hi < lo) — likely an algorithm bug")
    n == 1 && return lo
    k   = nbits(n)
    val = read_bits_int!(rng, k)
    val >= n && throw(OutOfRangeException())
    rng.coeff *= Q(2)^k // n
    lo + val
end

# Metropolis acceptance: mirrors  RandomReal[] < Piecewise[{{1, dE<=0}}, exp(-beta*exponent)].
# Always reads exactly one bit (always-read policy) so the symbolic decision
# tree matches Mathematica's; the symbolic factor is recorded, not evaluated.
# `exponent` defaults to dE (standard Metropolis); pass a scaled form for
# variants such as exp(-beta*dE/2).
function metropolis!(rng::BitSeqRNG, dE::LinForm; exponent::LinForm = dE)::Bool
    (any_tau_dep(dE) || any_tau_dep(exponent)) &&
        tau_violation!("acceptance threshold depends on an absolute (tau-dependent) position")
    cond = tau0_form(dE)
    expo = tau0_form(exponent)
    # dE identically 0  =>  clamp condition 0<=0 is True  =>  threshold 1  =>
    # always accept with no bit (matches acceptTestI's pVal>=hi shortcut).
    isempty(cond) && return true
    rng.pos += 1
    rng.pos > length(rng.bits) && throw(OutOfBitsException())
    accepted = rng.bits[rng.pos] == 1
    push!(rng.factors, MFactor(cond, expo, accepted))
    accepted
end

# ============================================================================
# SECTION 4 — BFS engine (per-state exhaustive path enumeration)
# ============================================================================
# A concrete state is a sorted vector of (row, col, type) integer triples.
const CState = Vector{NTuple{3,Int}}

intval(x::TauNum)::Int = (denominator(tau0(x)) == 1 ?
    Int(numerator(tau0(x))) : cant("non-integer position coordinate: $(tau0(x))"))

# Substitute tau = 0, apply PBC into 1..n, and sort -> canonical concrete state.
# tau=0 on a next-state position is always valid: a translation-invariant
# algorithm's next state is translation-COVARIANT (it shifts with the lattice),
# which is expected and must NOT be flagged.
function norm_state(s::PState, n::Int)::CState
    cs = NTuple{3,Int}[(mod(intval(p.r) - 1, n) + 1, mod(intval(p.c) - 1, n) + 1, p.t)
                       for p in s]
    sort!(cs)
    cs
end

struct Leaf
    next::CState
    coeff::Q
    factors::Vector{MFactor}
end

# Exhaustive BFS over all bit strings for one (tau-augmented) seed state.
# Returns the leaves; raises CantHandle if any path exceeds maxdepth (an
# incomplete tree would silently drop transitions and is treated as fatal).
function build_state_leaves(algo, seed::PState, n::Int, maxdepth::Int)::Vector{Leaf}
    leaves = Leaf[]
    queue  = Vector{Int}[Int[]]
    while !isempty(queue)
        bits = popfirst!(queue)
        rng  = BitSeqRNG(bits)
        try
            nxt = algo(rng, seed)::PState
            push!(leaves, Leaf(norm_state(nxt, n), rng.coeff, copy(rng.factors)))
        catch e
            if e isa OutOfBitsException
                if length(bits) < maxdepth
                    push!(queue, [bits; 0]); push!(queue, [bits; 1])
                else
                    cant("BFS incomplete: a random-number path reached maxdepth=$maxdepth " *
                         "before the algorithm returned a state; increase maxdepth")
                end
            elseif e isa OutOfRangeException
                # rejected (out-of-range) bit pattern in rejection sampling — drop
            else
                rethrow(e)
            end
        end
    end
    leaves
end

# ---- concrete evaluation of a leaf weight (sanity testing only) ----
# Threshold value of a factor at concrete couplings J (Float64) and beta.
function _eval_factor(f::MFactor, J::Dict{Atom,Float64}, beta::Float64)::Float64
    dE = sum(Float64(v) * get(J, a, 0.0) for (a, v) in f.cond; init = 0.0)
    g  = sum(Float64(v) * get(J, a, 0.0) for (a, v) in f.expo; init = 0.0)
    th = dE <= 0 ? 1.0 : exp(-beta * g)
    f.accepted ? th : (1.0 - th)
end
function eval_leaf(lf::Leaf, J::Dict{Atom,Float64}, beta::Float64)::Float64
    w = Float64(lf.coeff)
    for f in lf.factors; w *= _eval_factor(f, J, beta); end
    w
end

# ============================================================================
# SECTION 5 — State enumeration
# ============================================================================
# Ordered selections (k-permutations without repetition) of `items`.
function _kperms(items::Vector{T}, k::Int) where {T}
    out = Vector{T}[]
    n = length(items)
    used = falses(n)
    cur = Vector{T}(undef, k)
    function rec(d)
        if d > k
            push!(out, copy(cur)); return
        end
        for i in 1:n
            used[i] && continue
            used[i] = true; cur[d] = items[i]
            rec(d + 1)
            used[i] = false
        end
    end
    rec(1)
    out
end

# All distinct N-particle states for the given type multiset on the n x n torus.
function enumerate_states(types::Vector{Int}, n::Int)::Vector{CState}
    st  = sort(types)
    pos = [(r, c) for r in 1:n for c in 1:n]
    seen = Set{CState}()
    out  = CState[]
    for perm in _kperms(pos, length(st))
        cs = sort(NTuple{3,Int}[(perm[k][1], perm[k][2], st[k]) for k in 1:length(st)])
        if !(cs in seen)
            push!(seen, cs); push!(out, cs)
        end
    end
    out
end

# Theoretical combinatorial count: P(S,N) / prod(multiplicity!).
function theoretical_count(types::Vector{Int}, n::Int)
    S = n^2; N = length(types)
    num = prod(BigInt(S - k) for k in 0:N-1)
    den = prod(factorial(BigInt(c)) for c in values(_counts(types)))
    num ÷ den
end
function _counts(xs)
    d = Dict{Int,Int}()
    for x in xs; d[x] = get(d, x, 0) + 1; end
    d
end

# Concrete state -> tau-FREE PState (for energy evaluation of real states).
concrete_pstate(cs::CState)::PState = PState([Particle(TauNum(r), TauNum(c), t) for (r, c, t) in cs])
# Concrete state -> tau-AUGMENTED PState (BFS seed for the tau-check).
augmented_pstate(cs::CState)::PState = PState([aug_particle(r, c, t) for (r, c, t) in cs])

# ============================================================================
# SECTION 6 — T-matrix build (per-state BFS) + ergodicity reachability
# ============================================================================
# Runs the tau-augmented BFS from EVERY state (no orbit reduction), so the
# transition matrix is built directly and the tau-check covers all states.
struct BFSResult
    states        :: Vector{CState}
    idx           :: Dict{CState,Int}
    leaves        :: Vector{Vector{Leaf}}   # leaves[i] = leaves from state i
    tau_free      :: Bool
    tau_msg       :: String
end

function build_transitions(algo, energy, states::Vector{CState}, n::Int, maxdepth::Int)::BFSResult
    idx = Dict(cs => i for (i, cs) in enumerate(states))
    leaves = Vector{Vector{Leaf}}(undef, length(states))
    tau_free = true; tau_msg = ""
    for (i, cs) in enumerate(states)
        _TAU[] = false; _TAU_MSG[] = ""
        leaves[i] = build_state_leaves(algo, augmented_pstate(cs), n, maxdepth)
        if _TAU[]
            tau_free = false
            isempty(tau_msg) && (tau_msg = _TAU_MSG[])
        end
    end
    BFSResult(states, idx, leaves, tau_free, tau_msg)
end

# Reachability ergodicity: BFS over the directed transition graph from the seed.
function check_ergodicity(bfs::BFSResult, seed::CState)
    nS = length(bfs.states)
    adj = [Set{Int}() for _ in 1:nS]
    for i in 1:nS, lf in bfs.leaves[i]
        j = bfs.idx[lf.next]
        j != i && push!(adj[i], j)
    end
    seedi = bfs.idx[seed]
    seen = Set{Int}([seedi]); queue = [seedi]
    while !isempty(queue)
        u = popfirst!(queue)
        for v in adj[u]
            if !(v in seen); push!(seen, v); push!(queue, v); end
        end
    end
    (ergodic = length(seen) == nS, reached = length(seen), total = nS)
end

# ============================================================================
# SECTION 7 — Detailed-balance check (compact structure + Phase 2 + Phase 3)
# ============================================================================
# Mirrors dbc_core.wl's $dbcExportCompact + julia_poc/ebe.jl, but built fully
# in-memory and generalised to EXACT rationals (Rational{BigInt}) throughout,
# so half-integer Boltzmann exponents (e.g. exp(-beta*dE/2)) are handled exactly
# instead of triggering a Mathematica fallback.

using HiGHS

# Canonical hashable key for a RatForm (sorted (atom,coeff) pairs).
_form_key(f::RatForm) = Tuple(sort([(a, v) for (a, v) in f]))
# Canonical hashable key for a whole leaf weight (coeff + ordered factors).
_weight_key(lf::Leaf) =
    (lf.coeff, Tuple((_form_key(f.cond), _form_key(f.expo), f.accepted) for f in lf.factors))

vecof(f::RatForm, aidx::Dict{Atom,Int}, nA::Int) =
    (v = zeros(Q, nA); for (a, c) in f; v[aidx[a]] = c; end; v)

# ---- build all the index structures from the BFS result -------------------
struct DBModel
    atoms          :: Vector{Atom}
    aidx           :: Dict{Atom,Int}
    conds          :: Vector{RatForm}          # allConds (each is cond_lhs of cond_lhs<=0)
    cond_idx       :: Dict{Any,Int}
    cond_eff_lhs   :: Vector{Vector{Q}}        # sigma=1 iff eff_lhs . J >= 0
    cond_is_strict :: Vector{Bool}
    energy_coeffs  :: Vector{Vector{Q}}        # per state
    uweights       :: Vector{Leaf}             # representative leaf per unique weight
    uw_idx         :: Dict{Any,Int}
    uw_active      :: Vector{Vector{Int}}      # active cond indices per unique weight
    uw_cases       :: Vector{Dict{Vector{Int},Vector{Tuple{Q,Vector{Q}}}}}  # sigma(active)->terms
    pairs          :: Vector{Tuple{Int,Int}}   # (i<j) communicating pairs
    ij_srcs        :: Vector{Vector{Int}}      # unique-weight idx per i->j leaf
    ji_srcs        :: Vector{Vector{Int}}      # unique-weight idx per j->i leaf
end

function build_dbmodel(bfs::BFSResult, energy)::DBModel
    nS = length(bfs.states)

    # --- collect conditions (dedup by RatForm) ---
    conds = RatForm[]; cond_idx = Dict{Any,Int}()
    function cidx!(f::RatForm)
        k = _form_key(f)
        haskey(cond_idx, k) && return cond_idx[k]
        push!(conds, f); cond_idx[k] = length(conds); length(conds)
    end
    for i in 1:nS, lf in bfs.leaves[i], fac in lf.factors
        cidx!(fac.cond)
    end

    # --- atom list: union over conditions, exponents, and state energies ---
    state_energy = Vector{RatForm}(undef, nS)
    for i in 1:nS
        state_energy[i] = tau0_form(energy(concrete_pstate(bfs.states[i])))
    end
    atomset = Set{Atom}()
    for f in conds, a in keys(f); push!(atomset, a); end
    for i in 1:nS, lf in bfs.leaves[i], fac in lf.factors, a in keys(fac.expo); push!(atomset, a); end
    for i in 1:nS, a in keys(state_energy[i]); push!(atomset, a); end
    atoms = sort(collect(atomset))
    aidx  = Dict(a => i for (i, a) in enumerate(atoms)); nA = length(atoms)

    k = length(conds)
    cond_eff_lhs = Vector{Vector{Q}}(undef, k)
    cond_is_strict = fill(false, k)   # all conditions here are `<= 0` (LessEqual)
    for i in 1:k
        # LessEqual: cond TRUE iff cond_lhs <= 0 iff (-cond_lhs).J >= 0
        cond_eff_lhs[i] = -vecof(conds[i], aidx, nA)
    end

    energy_coeffs = [vecof(state_energy[i], aidx, nA) for i in 1:nS]

    # --- unique weights + compact term tables ---
    uweights = Leaf[]; uw_idx = Dict{Any,Int}()
    function widx!(lf::Leaf)
        key = _weight_key(lf)
        haskey(uw_idx, key) && return uw_idx[key]
        push!(uweights, lf); uw_idx[key] = length(uweights); length(uweights)
    end

    # --- ordered transitions ---
    trans = Dict{Tuple{Int,Int},Vector{Int}}()   # (i,j)->[uw idx per leaf]
    for i in 1:nS, lf in bfs.leaves[i]
        j = bfs.idx[lf.next]
        j == i && continue
        push!(get!(trans, (i, j), Int[]), widx!(lf))
    end

    uw_active, uw_cases = _compact_weights(uweights, cond_idx, aidx, nA, k)

    # --- communicating pairs ---
    pairset = Set{Tuple{Int,Int}}()
    for (i, j) in keys(trans)
        push!(pairset, (min(i, j), max(i, j)))
    end
    pairs = sort(collect(pairset))
    ij_srcs = [get(trans, (a, b), Int[]) for (a, b) in pairs]
    ji_srcs = [get(trans, (b, a), Int[]) for (a, b) in pairs]

    DBModel(atoms, aidx, conds, cond_idx, cond_eff_lhs, cond_is_strict,
            energy_coeffs, uweights, uw_idx, uw_active, uw_cases,
            pairs, ij_srcs, ji_srcs)
end

# Build, per unique weight, the active condition indices and the term table
# for every truth assignment of those conditions (sigma-substitution).
function _compact_weights(uweights, cond_idx, aidx, nA, k)
    uw_active = Vector{Vector{Int}}(undef, length(uweights))
    uw_cases  = Vector{Dict{Vector{Int},Vector{Tuple{Q,Vector{Q}}}}}(undef, length(uweights))
    for (wi, lf) in enumerate(uweights)
        active = sort(unique(Int[cond_idx[_form_key(f.cond)] for f in lf.factors]))
        uw_active[wi] = active
        ka = length(active)
        apos = Dict(c => p for (p, c) in enumerate(active))
        cases = Dict{Vector{Int},Vector{Tuple{Q,Vector{Q}}}}()
        for mask in 0:(2^ka - 1)
            asg = [(mask >> (p - 1)) & 1 for p in 1:ka]   # value per active cond
            terms = Tuple{Q,Vector{Q}}[(lf.coeff, zeros(Q, nA))]
            for f in lf.factors
                ci = cond_idx[_form_key(f.cond)]
                sig = asg[apos[ci]]                       # 1 = cond TRUE
                g = vecof(f.expo, aidx, nA)
                if f.accepted
                    if sig == 0                            # threshold exp(-beta*g)
                        terms = [(c, e .+ g) for (c, e) in terms]
                    end                                    # sig==1 -> factor 1
                else
                    if sig == 1                            # 1 - 1 = 0
                        terms = Tuple{Q,Vector{Q}}[]
                    else                                   # 1 - exp(-beta*g)
                        terms = vcat(terms, [(-c, e .+ g) for (c, e) in terms])
                    end
                end
            end
            cases[asg] = _combine_terms(terms, nA)
        end
        uw_cases[wi] = cases
    end
    uw_active, uw_cases
end

function _combine_terms(terms, nA)
    acc = Dict{Vector{Q},Q}()
    for (c, e) in terms
        acc[e] = get(acc, e, Q(0)) + c
    end
    [(c, e) for (e, c) in acc if c != 0]
end

# ---- Phase 2: chamber enumeration via HiGHS LP BFS + degenerate filter ----
function _is_feasible(sigma::Vector{Int}, eff::Vector{Vector{Q}},
                     strict::Vector{Bool}, nA::Int, eps::Float64)::Bool
    isempty(sigma) && return true
    h = HiGHS.Highs_create()
    HiGHS.Highs_setBoolOptionValue(h, "output_flag", false)
    for _ in 1:nA
        HiGHS.Highs_addVar(h, -1e10, 1e10)
    end
    inds = Int32[j - 1 for j in 1:nA]
    for i in eachindex(sigma)
        s = sigma[i] == 1 ? 1.0 : -1.0
        coeffs = Float64[s * Float64(eff[i][j]) for j in 1:nA]
        lb = ((sigma[i] == 1) == strict[i]) ? eps : 0.0
        HiGHS.Highs_addRow(h, lb, 1e30, nA, inds, coeffs)
    end
    HiGHS.Highs_run(h)
    st = HiGHS.Highs_getModelStatus(h)
    HiGHS.Highs_destroy(h)
    st == 7   # kOptimal
end

function _filter_degenerate(feasible, eff, strict, k)
    contra_ff = Tuple{Int,Int}[]; contra_tt = Tuple{Int,Int}[]
    for i in 1:k, j in (i+1):k
        if eff[i] == -eff[j]
            if strict[i] && strict[j]
                push!(contra_ff, (i, j))
            elseif !strict[i] && !strict[j]
                push!(contra_tt, (i, j))
            end
        end
    end
    (isempty(contra_ff) && isempty(contra_tt)) && return feasible
    filter(feasible) do s
        !any(s[i] == 0 && s[j] == 0 for (i, j) in contra_ff) &&
        !any(s[i] == 1 && s[j] == 1 for (i, j) in contra_tt)
    end
end

function enumerate_chambers(m::DBModel)::Vector{Vector{Int}}
    k = length(m.conds); nA = length(m.atoms)
    k == 0 && return [Int[]]
    # generic interior witness: J*_a = 100^a guarantees no real condition is 0.
    Jstar = Q[Q(100)^a for a in 1:nA]
    initial = Int[ (sum(m.cond_eff_lhs[i][a] * Jstar[a] for a in 1:nA) >= 0) ? 1 : 0
                   for i in 1:k ]
    eps = 1e-6
    visited = Set{Vector{Int}}([copy(initial)])
    feasible = [copy(initial)]; queue = [copy(initial)]
    while !isempty(queue)
        s = popfirst!(queue)
        for i in 1:k
            s2 = copy(s); s2[i] = 1 - s2[i]
            s2 in visited && continue
            push!(visited, copy(s2))
            if _is_feasible(s2, m.cond_eff_lhs, m.cond_is_strict, nA, eps)
                push!(feasible, copy(s2)); push!(queue, copy(s2))
            end
        end
    end
    _filter_degenerate(feasible, m.cond_eff_lhs, m.cond_is_strict, k)
end

# ---- Phase 3: exact rational DB check over all chambers --------------------
# Project a full chamber sigma onto a unique weight's active conditions.
function _terms_for(m::DBModel, wi::Int, sigma::Vector{Int})
    asg = Int[sigma[c] for c in m.uw_active[wi]]
    m.uw_cases[wi][asg]
end

# Returns (pass, violations) where each violation is (i, j, chamber_index).
function run_db_check(m::DBModel)
    chambers = enumerate_chambers(m)
    violations = Tuple{Int,Int,Int}[]
    for (ridx, sigma) in enumerate(chambers)
        for p in eachindex(m.pairs)
            (isempty(m.ij_srcs[p]) && isempty(m.ji_srcs[p])) && continue
            i, j = m.pairs[p]
            ei = m.energy_coeffs[i]; ej = m.energy_coeffs[j]
            groups = Dict{Vector{Q},Q}()
            for wi in m.ij_srcs[p], (c, g) in _terms_for(m, wi, sigma)
                key = g .+ ei
                groups[key] = get(groups, key, Q(0)) + c
            end
            for wi in m.ji_srcs[p], (c, g) in _terms_for(m, wi, sigma)
                key = g .+ ej
                groups[key] = get(groups, key, Q(0)) - c
            end
            if any(c != 0 for c in values(groups))
                push!(violations, (i, j, ridx))
            end
        end
    end
    (isempty(violations), violations, length(chambers))
end
