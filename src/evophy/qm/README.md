# evophy.qm — Schrödinger / quantum branch

The `qm-schrodinger` branch extends evophy's GP/analytical-discovery approach
from classical Hamiltonian mechanics (2D Kepler) to the **1D time-independent
Schrödinger equation**

```
  -(ħ²/2m) ψ''(x) + V(x) ψ(x) = E ψ(x)
```

The long-term target is **nuclear-reaction-relevant potentials** with no
closed-form spectrum. Following the pattern that worked for Kepler, we start
by validating the ground-truth oracle on **solvable benchmarks** before
pointing search at open problems.

## What's here so far

- **`schrodinger.clj`** — the ground-truth oracle (analogue of the symplectic
  integrator for Kepler). A pure-Clojure **Numerov + shooting** eigenvalue
  solver: for a given `V(x)` it returns bound-state energies `E_n` and unit-
  normalized wavefunctions `ψ_n(x)`.
  - Eigenvalues are bracketed by counting nodes in the *classically-allowed*
    region `V ≤ E` (robust against runaway tails), then refined on the sign
    of the far-boundary value `ψ(end)` — the classic shooting criterion,
    accurate for both hard walls and exponentially-decaying states.
  - Benchmarks with known spectra: `infinite-well` (`E_n = n²π²ħ²/2mL²`),
    `harmonic` (`E_n = ħω(n+½)`), and `finite-well` (transcendental closed form
    via `finite-well-levels`).
  - **`woods-saxon`** — the first genuine target: a smooth nuclear mean-field
    well `V(x) = -V₀/(1+e^((|x|-R)/a))` with **no closed-form spectrum**,
    validated only by grid convergence and internal consistency (ordering,
    node counts), not against an analytic answer.

- **`fitness.clj`** — the scoring layer (roadmap item 2). Grades what a search
  proposes against the oracle, keeping the classical branch's "approximately
  right, and cheap counts too" spirit rather than exact-or-nothing:
  - **Energy / spectrum scoring** — `energy-score`, `spectrum-score` give
    partial credit via `1/(1 + (rel/tol)²)` (the same reward shape the
    classical branch uses), rewarding both the values and the level count.
  - **Wavefunction fidelity** — `overlap` / `wavefn-score` compute the
    phase-independent overlap `|⟨ψ_trial|ψ_ref⟩|` against oracle eigenstates.
  - **Variational objective (oracle-free)** — `rayleigh-quotient` evaluates
    `⟨ψ|H|ψ⟩/⟨ψ|ψ⟩`, a rigorous *upper bound* on the ground-state energy, so a
    trial ψ can be scored on physics alone; `project-out` orthogonalizes
    against lower states to bound excited levels.
  - **Accuracy vs. cost** — `pareto-dominates?` / `pareto-front` keep the two
    objectives separate; `scalar-fitness` collapses them only when a search
    wants one number. `score-strategy` runs a whole procedure and scores its
    spectrum + reported work end-to-end.

- **`search.clj`** — the search layer (roadmap item 3). A fresh, small GP
  loop over trial-ψ expression trees, reusing only the classical branch's
  *pattern* (truncation elitism + mutation/crossover + immigrants), not any
  of its Kepler code:
  - **Genome** — a single s-expression over `x` (plus the potential `V`
    itself, passed in as a live argument so a compiled fn is valid for any
    benchmark), with a small generic op set (`+ - * qdiv`, `qexp qneg
    qsquare qabs qsqrt qsin qcos`, `V`) — no orbital primitives.
  - **Fitness is oracle-free** — scores against `fitness.clj`'s
    `variational-fitness` (Rayleigh quotient), a rigorous upper bound on the
    true energy, so search climbs toward truth using only the potential and
    the variational principle, never the answer. Excited states come from
    sequential Gram-Schmidt deflation (`project-out`) against the levels
    already found.
  - `evolve-level` searches one level; `evolve-spectrum` chains levels
    bottom-up. On the harmonic oscillator this rediscovers the exact ground
    and first-excited states (`exp(-V(x))`, `x*exp(-V(x))`) from scratch.

Run the validation and the search/fitness demos:

```bash
lein run -m evophy.qm.schrodinger    # oracle benchmarks
lein run -m evophy.qm.fitness        # scoring demo on the harmonic oscillator
lein run -m evophy.qm.search         # GP search discovers E_0..E_2 from scratch
lein test evophy.qm.schrodinger-test evophy.qm.fitness-test evophy.qm.search-test
```

## Roadmap

1. ~~**More benchmarks** — finite square well, then a nuclear-like potential
   (Woods–Saxon).~~ **Done.** Solver validated on infinite well (~1e-11),
   harmonic (~1e-10), finite well (transcendental, grid-limited by the
   potential's discontinuity), and Woods–Saxon (grid-converged, no closed form).
2. ~~**Fitness layer** — score trial energies / wavefunctions / *strategies*
   against the oracle. Keep the "approximately right and cheap counts too"
   spirit (partial credit, Pareto accuracy-vs-cost), not exact-or-nothing.~~
   **Done.** See `fitness.clj`: partial-credit energy/spectrum scores,
   oracle-free variational objective (Rayleigh quotient + orthogonal
   projection), wavefunction fidelity, and Pareto accuracy-vs-cost.
3. ~~**Search layer** — GP over trial-ψ expression trees and/or procedural
   refinement strategies (Newton on the boundary residual, variational
   trial functions, shooting-step rules). Reuse philosophy from the classical
   branch, not its Kepler-specific primitives.~~ **Done** (first pass). See
   `search.clj`: variational (oracle-free) GP over trial-ψ trees, ground
   state + sequential deflation for excited states. Not yet explored:
   procedural refinement strategies, and pointing the search at Woods-Saxon
   (no closed form to validate against, but nothing stops it from running).

Nothing from the classical `evophy.core` pipeline (phase-space representation,
orbital primitives, coordinate charts) transfers directly — this is a fresh
pipeline sharing only the method.
