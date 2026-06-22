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
| `--de-driven` | Fitness from known equations of motion (ODE residual / orbit invariance), not trajectory matching |
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
| `--mcts-simulations N` | MCTS rollouts when using `evophy.mcts` (default 64) |
| `--mcts-inject N` | (CLI only) |

Resume training: run `lein run` again without `--fresh`.

## Individuals and strategies

Each individual holds one or more **laws** in `:laws`. Composite fitness is the **worst** law score (all must be good).

| Kind | Role | Fitness (data-driven) |
|------|------|------------------------|
| `:differential` | Rate laws \(\dot q_x, \dot q_y, \dot p_x, \dot p_y\) from current `(qx, qy, px, py, m, alpha)` | Symplectic one-step / rollout error |
| `:analytical` | Trajectories \(q(t), p(t)\) from `(t, ICs, m, alpha)` | Short-horizon forecast error |
| `:conserved` | Scalar invariant \(C(q,p,m,\alpha)\) | Low CoV along orbit **and** \(\nabla C \cdot f \approx 0\) |

With `--de-driven`, analytical laws score ODE residual on phase samples; conserved laws score **temporal constancy along integrated orbits** (they are invariants, not solutions of the motion DE). Differential laws score 0 (the DE is already known).

Each scenario supplies its own `m` and `alpha` when evaluating fitness.

## Coordinate charts

Integration stays Cartesian; evolved laws may use the chart best suited to each scenario's ICs.

- **Cartesian** — `qx, qy, px, py` (and ICs `q0x, q0y, p0x, p0y`).
- **Polar** — `r, theta, pr, ptheta` (and ICs `r0, theta0, pr0, ptheta0`).

Each scenario dataset carries a `:chart` (polar when the IC is axisymmetric in the \(x\)-axis plane, e.g. circle, heavy-m, strong-g). **During fitness, a law uses that scenario chart unless the law sets `:chart` explicitly** — so circle runs evolve in polar without tagging every law. New immigrants default to the dominant chart of the training batch.

Laws with the wrong expression keys for the effective chart score 0 (e.g. Cartesian-only analytical laws on a polar scenario).

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
