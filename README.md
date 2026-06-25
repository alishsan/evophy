# evophy

Evolve symbolic **2D gravitational** dynamics from synthetic trajectories. Ground truth uses a symplectic integrator for \(H = (p_x^2 + p_y^2)/(2m) - \alpha/r\). Search is **genetic programming** over composable expression blocks, with optional MCTS helpers for analytical genomes (`evophy.mcts`).

Repository: [https://github.com/alishsan/evophy](https://github.com/alishsan/evophy)

## Requirements

- [Leiningen](https://leiningen.org/) 2.x
- JDK 11+

## Quick start

```bash
lein test
lein run
```

`lein run` trains on several 2D scenarios (different ICs, \(m\), \(\alpha\)), saves the population to `data/population.edn`, and prints the top distinct genomes plus per-scenario MSE.

Example run focused on conserved invariants with physics-driven fitness (analytical ODE residual + orbit invariance):

```bash
lein run -- --fresh --de-driven --strategy conserved --generations 50 --population-size 50 --scenario-seed 42
```

Use `--fresh` after checkpoint format changes or when you want a clean population.

### CLI

| Flag | Meaning |
|------|---------|
| `--fresh` | Ignore checkpoint; new random population |
| `--seed` | Inject physics seed individuals into the first generation |
| `--de-driven` | Fitness from known physics: analytical short-horizon orbit fit + conserved invariance along orbits (not data-driven scenario matching) |
| `--strategy NAME` | Restrict immigrants to `analytical`, `differential`, or `conserved` (default: mixed) |
| `--generations N` | Generations this run (default 50) |
| `--population-size N` | Population size (default 50) |
| `--population PATH` | Checkpoint file (default `data/population.edn`) |
| `--fixed-scenarios` | Train on the fixed reference scenario set each generation |
| `--random-scenarios` | Sample random scenarios each generation (default) |
| `--scenario-samples N` | Scenarios per fitness batch when random (default 32) |
| `--scenario-seed N` | RNG seed for scenario sampling |
| `--fitness-aggregate min\|percentile` | Aggregate fitness across scenarios (default `min`) |
| `--fitness-percentile P` | Percentile when aggregate is `percentile` (default 10) |
| `--no-guess` | Disable guess-style subtree mutations |
| `--prompt-each-generation` | Pause between generations (Enter / q) |
| `--no-mcts` | (CLI only; main loop is GP) |
| `--no-mcts-repair` | Disable adversarial MCTS slot repair immigrants |
| `--mcts-repair-simulations N` | MCTS sims per repair (default 48) |
| `--mcts-repair-inject N` | Max repair immigrants per generation (default 1; more during stagnation) |
| `--no-mcts-mutate` | Disable MCTS slot repair during analytical mutation |
| `--mcts-mutate-rate R` | Fraction of analytical mutates that try MCTS repair (default 0.2) |
| `--mcts-mutate-simulations N` | MCTS sims per mutation attempt (default 48) |
| `--mcts-simulations N` | MCTS rollouts when using `evophy.mcts` (default 64) |
| `--mcts-inject N` | (CLI only) |

Resume training: run `lein run` again without `--fresh`.

## Individuals and strategies

Each individual holds one or more **laws** in `:laws`. Composite fitness is the **worst** law score (all must be good).

| Kind | Role | Fitness (data-driven) |
|------|------|------------------------|
| `:differential` | Rate laws \(\dot q_x, \dot q_y, \dot p_x, \dot p_y\) from current `(qx, qy, px, py, m, alpha)` | Symplectic one-step / rollout error |
| `:analytical` | Trajectories \(q(t), p(t)\) from `(t, ICs, m, alpha)` | Short-horizon forecast error (matching **domain** or `e/if` branch) |
| `:conserved` | Scalar invariant \(C(q,p,m,\alpha)\) | Low CoV along orbit **and** \(\nabla C \cdot f \approx 0\) |

Analytical validity can be declared two ways:

1. **`:domain` tag** — `:bound` (\(E<0\)), `:unbound` (\(E>0\)), or `:any` (legacy). Fitness and MSE run only on scenarios in that regime; mismatches are **n/a**, not penalized. Use `--domain-filter` with `--de-driven --strategy analytical` to keep this mode.
2. **`(e/if (neg? energy) bound-branch unbound-branch)`** inside expressions — scored on every scenario; `energy` is \(H(q_0,p_0)\) pre-bound at compile time.

**Both regimes (default for `--de-driven --strategy analytical`):** fitness uses **all** reference scenarios; `:domain` tags are ignored; every analytical slot must be `(e/if (neg? energy) bound-arm unbound-arm)` with no bogus `e/if` tests. **Per-regime arm fitness:** bound arms are scored only on bound scenarios, unbound arms on unbound; overall fitness is `min(bound, unbound)`. Seeds are wrapped automatically. Opt out with `--domain-filter`.

With `--de-driven`, analytical laws score short-horizon trajectory fit on integrated reference orbits; conserved laws score **temporal constancy along integrated orbits** (they are invariants, not solutions of the motion DE). Differential laws score 0 (the DE is already known).

Each scenario supplies its own `m` and `alpha` when evaluating fitness.

## Coordinate charts

Integration stays Cartesian; evolved laws may use the chart best suited to each scenario's ICs.

- **Cartesian** — `qx, qy, px, py` (and ICs `q0x, q0y, p0x, p0y`).
- **Polar** — `r, theta, pr, ptheta` (and ICs `r0, theta0, pr0, ptheta0`).

Each scenario dataset carries a `:chart` (polar when the IC is axisymmetric in the \(x\)-axis plane, e.g. circle, heavy-m, strong-g). **Analytical laws are evaluated in the chart implied by their expression keys** (`qx`/`qy`… or `r`/`theta`…), using enriched trajectory data — so a Cartesian law still scores on a nominally polar scenario. Conserved laws still use the scenario chart unless the law sets `:chart` explicitly.

## Tests

```bash
lein test
```

## Uberjar (optional)

```bash
lein uberjar
java -jar target/uberjar/evophy-*-standalone.jar
```

## License

**Eclipse Public License 2.0** — see [`LICENSE`](LICENSE).
