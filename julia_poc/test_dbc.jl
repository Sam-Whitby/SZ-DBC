# ============================================================================
# test_dbc.jl  —  Regression suite for the pure-Julia SZ-DBC checker
# ============================================================================
# Proves the individual pieces (TauNum, bit-exact selection weights, symbolic
# leaf weights, tau detection) and runs two fast end-to-end checks (one DB PASS,
# one DB FAIL) so a false positive OR a false negative would fail the suite.
#
# Run:  julia --project=julia_poc julia_poc/test_dbc.jl
# ============================================================================

include(joinpath(@__DIR__, "dbc.jl"))
using Test, Random

@testset "SZ-DBC Julia checker" begin

@testset "TauNum: linear tracking + nonlinear taint" begin
    a = tau_r_aug(3); b = tau_r_aug(1)
    @test is_tau_free(a - b) && tau0(a - b) == 2      # difference cancels tau
    @test !is_tau_free(a^2) && tau0(a^2) == 9          # square taints
    @test !is_tau_free(a + 5) && tau0(a + 5) == 8      # affine keeps tau
    e = a * 2
    @test !is_tau_free(e) && e.cr == 2 && tau0(e) == 6
    @test is_tau_free((a - b) * tau_const(4))          # free*free stays free
    @test !is_tau_free(tau_r_aug(1) * tau_c_aug(1))    # tau_r*tau_c taints
end

@testset "nbits matches IntegerLength[n-1,2]" begin
    @test [nbits(n) for n in 1:9] == [0,1,2,2,3,3,3,3,4]
end

@testset "rejection-sampling weights are exact (1/n), n leaves" begin
    function choice_weights(n)
        total = Q(0); leaves = 0; queue = [Int[]]
        while !isempty(queue)
            bits = popfirst!(queue); rng = BitSeqRNG(bits)
            try
                rand_choice_index!(rng, n); total += rng.coeff; leaves += 1
            catch e
                e isa OutOfBitsException ? (push!(queue,[bits;0]); push!(queue,[bits;1])) :
                e isa OutOfRangeException ? nothing : rethrow(e)
            end
        end
        (leaves, total)
    end
    for n in 1:9
        l, w = choice_weights(n)
        @test l == n && w == 1
    end
end

# Shared mini single-particle Metropolis on an n x n torus (pairwise energy).
DISPS8 = [(dx,dy) for dx in -1:1 for dy in -1:1 if (dx,dy)!=(0,0)]
function mk_energy(n, maxd2)
    (state::PState) -> begin
        lf = LinForm()
        for i in 1:length(state), j in (i+1):length(state)
            d2 = pbc_d2(state[i], state[j], n); (0 < d2 <= maxd2) || continue
            addcoef!(lf, Jc(state[i].t, state[j].t, d2), TauNum(1))
        end
        lf
    end
end
function mk_metropolis(n, maxd2, energy)
    (rng, state::PState) -> begin
        pidx = rand_choice_index!(rng, length(state)); p = state[pidx]
        (dr,dc) = rand_choice!(rng, DISPS8); newp = Particle(p.r+dr, p.c+dc, p.t)
        rest = state[setdiff(1:length(state), pidx)]
        for q in rest; same_site(q, newp, n) && return state; end
        ns = vcat(rest, [newp]); dE = linsub(energy(ns), energy(state))
        metropolis!(rng, dE) ? ns : state
    end
end

@testset "leaf weights sum to 1 at random coupling points" begin
    n=3; energy = mk_energy(n,2); step = mk_metropolis(n,2,energy)
    seed = augmented_pstate(sort([(1,1,1),(2,2,2),(3,3,3)]))   # spread (uphill moves)
    _TAU[] = false
    leaves = build_state_leaves(step, seed, n, 22)
    @test _TAU[] == false
    Random.seed!(1)
    atoms = unique(vcat([collect(keys(f.cond)) for l in leaves for f in l.factors]...))
    for _ in 1:5
        J = Dict{Atom,Float64}(a => rand()*2-1 for a in atoms)
        @test abs(sum(eval_leaf(l, J, 1.0) for l in leaves) - 1.0) < 1e-9
    end
end

@testset "tau detection: no false positive, no false negative" begin
    n=2
    # absolute-position (quadratic field) energy MUST flag
    energy_qf = function(state::PState)
        lf = mk_energy(n,2)(state)
        fld = TauNum(0); for p in state; fld = fld + p.r^2; end
        addcoef!(lf, Xparam(:fieldH), fld); lf
    end
    qf = mk_metropolis(n,2,energy_qf)
    seed = augmented_pstate(sort([(1,1,1),(1,2,2)]))
    _TAU[]=false; build_state_leaves(qf, seed, n, 22)
    @test _TAU[] == true
    # pairwise-only energy MUST NOT flag
    energy_p = mk_energy(n,2); pstep = mk_metropolis(n,2,energy_p)
    _TAU[]=false; build_state_leaves(pstep, seed, n, 22)
    @test _TAU[] == false
end

@testset "end-to-end DB PASS (correct metropolis, nGrid=2)" begin
    n=2; energy = mk_energy(n,2); step = mk_metropolis(n,2,energy)
    states = enumerate_states([1,2], n)
    bfs = build_transitions(step, energy, states, n, 22)
    pass, viol, nch = run_db_check(build_dbmodel(bfs, energy))
    @test bfs.tau_free
    @test pass                       # no false NEGATIVE on a correct algorithm
end

@testset "end-to-end DB FAIL (variable-pool 4-hop)" begin
    n=3; const_bvp = [(-1,0),(1,0),(0,-1),(0,1)]
    energy0 = (state::PState) -> LinForm()
    bvp = function(rng, state::PState)
        pidx = rand_choice_index!(rng, length(state)); p = state[pidx]
        rest = state[setdiff(1:length(state), pidx)]
        valid = [(d1,d2) for (d1,d2) in const_bvp
                 if !any(same_site(q, Particle(p.r+d1,p.c+d2,p.t), n) for q in rest)]
        isempty(valid) && return state
        (d1,d2) = rand_choice!(rng, valid)
        vcat(rest, [Particle(p.r+d1, p.c+d2, p.t)])
    end
    states = enumerate_states([1,2], n)
    bfs = build_transitions(bvp, energy0, states, n, 22)
    pass, viol, nch = run_db_check(build_dbmodel(bfs, energy0))
    @test bfs.tau_free
    @test !pass                      # no false POSITIVE on a broken algorithm
end

@testset "exact rationals: half-beta exponent is detected as FAIL" begin
    n=3; energy = mk_energy(n,2)
    halfstep = function(rng, state::PState)
        pidx = rand_choice_index!(rng, length(state)); p = state[pidx]
        (dr,dc) = rand_choice!(rng, DISPS8); newp = Particle(p.r+dr, p.c+dc, p.t)
        rest = state[setdiff(1:length(state), pidx)]
        for q in rest; same_site(q, newp, n) && return state; end
        ns = vcat(rest, [newp]); dE = linsub(energy(ns), energy(state))
        half = LinForm(a => c*(1//2) for (a,c) in dE)
        metropolis!(rng, dE; exponent=half) ? ns : state
    end
    states = enumerate_states([1,2,3], n)
    bfs = build_transitions(halfstep, energy, states, n, 22)
    pass, viol, nch = run_db_check(build_dbmodel(bfs, energy))
    @test !pass
end

end
