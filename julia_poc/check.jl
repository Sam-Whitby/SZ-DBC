#!/usr/bin/env julia
# ============================================================================
# check.jl  —  Pure-Julia SZ-DBC checker CLI (no Mathematica)
# ============================================================================
# Usage:
#     julia --project=julia_poc julia_poc/check.jl j_examples/<algorithm>.jl [options]
#
# Options:
#     -maxdepth N     Maximum BFS bit depth (default 22)
#
# The algorithm file (a Julia translation of an examples/*.wl algorithm) must
# define, as plain top-level bindings:
#
#     const NGRID            ::Int
#     const MAXD2            ::Int            # max squared distance for energy terms
#     const PARTICLE_TYPES   ::Vector{Int}
#     const SYMMETRY_GROUP   ::Vector{String}   # e.g. ["translation"]
#     energy(state::PState)  ::LinForm        # symbolic energy (linear in atoms)
#     algorithm(rng, state::PState) ::PState  # one MCMC step, using rng primitives
#
# The checker runs the ENTIRE pipeline in Julia:
#     enumerate states -> tau-augmented BFS (translational-invariance check +
#     exact path enumeration) -> ergodicity -> detailed balance (HiGHS chamber
#     enumeration + exact rational grouping).
#
# It returns exit code 0 iff: state count OK, ergodic, DB satisfied, and
# (translation invariance holds OR translation was not declared). Any
# unsupported construct raises a hard error and aborts (never a wrong verdict).
# ============================================================================

include(joinpath(@__DIR__, "dbc.jl"))

const SEP  = "="^64
const SEP2 = "-"^64

function parse_args(argv)
    isempty(argv) && (println("Usage: julia check.jl <algorithm.jl> [-maxdepth N]"); exit(1))
    algfile = argv[1]
    isfile(algfile) || (println("ERROR: file not found: ", algfile); exit(1))
    maxdepth = 22
    i = 2
    while i <= length(argv)
        if argv[i] == "-maxdepth"
            maxdepth = parse(Int, argv[i+1]); i += 2
        else
            println("ERROR: unknown option ", argv[i]); exit(1)
        end
    end
    (algfile, maxdepth)
end

# Canonical seed state (matches the .wl $seedState: first N row-major sites,
# sorted types).  Used as the ergodicity-reachability start.
function canonical_seed(types::Vector{Int}, n::Int)::CState
    st  = sort(types)
    pos = [(r, c) for r in 1:n for c in 1:n]
    sort(NTuple{3,Int}[(pos[k][1], pos[k][2], st[k]) for k in 1:length(st)])
end

function run_checker(algfile, maxdepth, n, types, symgrp, algo, energy)
    "D4" in symgrp && (println("ERROR: D4 symmetry is not supported by the Julia checker ",
                               "(translation-only). Aborting rather than ignoring it."); exit(1))
    want_trans = "translation" in symgrp

    println(SEP)
    println("  SZ-DBC (Julia)  —  ", basename(algfile))
    println(SEP)
    println("  nGrid      : ", n)
    println("  particles  : ", types)
    println("  symmetry   : ", symgrp)

    seed = canonical_seed(types, n)
    println("  seed state : ", seed)
    println()

    # ---- Step 1: state enumeration + combinatorial count ----
    println(SEP2); println("  Step 1: State enumeration"); println(SEP2)
    t1 = @elapsed states = enumerate_states(types, n)
    theo = theoretical_count(types, n)
    count_ok = length(states) == theo
    @printf("  States found : %d  (%.2fs)\n", length(states), t1)
    println("  State count  : ", count_ok ? "OK" : "MISMATCH",
            "  (found ", length(states), ", theoretical ", theo, ")")

    # ---- Step 2: tau-augmented BFS (translational invariance + path enumeration) ----
    println(SEP2); println("  Step 2: tau-BFS (translational invariance + path enumeration)"); println(SEP2)
    local bfs
    try
        t2 = @elapsed (bfs = build_transitions(algo, energy, states, n, maxdepth))
        @printf("  BFS over %d states  (%.2fs)\n", length(states), t2)
    catch e
        if e isa CantHandle
            println("  ERROR: ", e.msg); exit(1)
        else
            rethrow(e)
        end
    end
    if bfs.tau_free
        println("  Translational : PASS  — tau cancels in all leaf weights")
    else
        println("  Translational : FAIL  — ", bfs.tau_msg)
        println("  (DB is still checked directly per state; translation reduction not assumed.)")
    end

    # ---- Step 3: ergodicity (reachability) ----
    println(SEP2); println("  Step 3: Ergodicity (reachability from seed)"); println(SEP2)
    t3 = @elapsed erg = check_ergodicity(bfs, seed)
    @printf("  Reachable : %d/%d  (%.2fs)\n", erg.reached, erg.total, t3)
    println("  Ergodicity : ", erg.ergodic ? "PASS" : "FAIL")

    # ---- Step 4: detailed balance ----
    println(SEP2); println("  Step 4: Detailed balance (HiGHS chambers + exact rational check)"); println(SEP2)
    local pass, viol, nch
    try
        t4 = @elapsed ((pass, viol, nch) = (run_db_check(build_dbmodel(bfs, energy))))
        @printf("  Chambers : %d   Time: %.2fs\n", nch, t4)
    catch e
        if e isa CantHandle
            println("  ERROR: ", e.msg); exit(1)
        else
            rethrow(e)
        end
    end
    if pass
        println("  Detailed bal. : PASS  — satisfied for all ", length(states), " states")
    else
        println("  Detailed bal. : FAIL  — ", length(viol), " violating (pair, chamber) record(s):")
        for v in first(viol, min(5, length(viol)))
            println("    s=", states[v[1]], "  t=", states[v[2]], "  (chamber ", v[3], ")")
        end
    end

    # ---- Summary ----
    println(SEP); println("  SUMMARY"); println(SEP)
    println("  Translational : ", want_trans ? (bfs.tau_free ? "PASS" : "FAIL") : "not declared")
    println("  State count   : ", count_ok ? "OK" : "MISMATCH")
    println("  Ergodicity    : ", erg.ergodic ? "PASS" : "FAIL")
    println("  Detailed bal. : ", pass ? "PASS" : "FAIL")
    println(SEP)

    all_ok = pass && erg.ergodic && count_ok && (bfs.tau_free || !want_trans)
    exit(all_ok ? 0 : 1)
end

# Parse args and load the translated algorithm at TOP LEVEL (before the pipeline
# runs) so the user's `algorithm`/`energy` methods are visible to the call site
# under Julia 1.12's stricter world-age rules.
const _ALGFILE, _MAXDEPTH = parse_args(ARGS)
include(abspath(_ALGFILE))
run_checker(_ALGFILE, _MAXDEPTH,
            Main.NGRID, Main.PARTICLE_TYPES, Main.SYMMETRY_GROUP,
            Main.algorithm, Main.energy)
