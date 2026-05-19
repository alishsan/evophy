# evophy

Clojure experiment: evolve symbolic **2D gravitational** dynamics from synthetic trajectories. Ground truth uses a symplectic integrator for \(H = (p_x^2 + p_y^2)/(2m) - \alpha/r\). Search combines **genetic programming** with optional **MCTS** (analytical genomes only).

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

### CLI

| Flag | Meaning |
|------|---------|
| `--fresh` | Ignore checkpoint; new random population |
| `--no-mcts` | GP only (no MCTS injections) |
| `--generations N` | Generations this run (default 50) |
| `--population-size N` | Population size (default 50) |
| `--population PATH` | Checkpoint file (default `data/population.edn`) |
| `--mcts-simulations N` | MCTS rollouts per injection (default 64) |
| `--mcts-inject N` | Analytical individuals added per generation (default 5) |

Resume training: run `lein run` again without `--fresh`.

## Genome strategies

- **`:differential`** — four rate laws \(\dot q_x, \dot q_y, \dot p_x, \dot p_y\) as functions of **current** `(qx, qy, px, py)` and parameters **`m`**, **`alpha`**; fitness = one-step Euler forecast error along each trajectory.
- **`:analytical`** — four laws \(q_x(t,\ldots), \ldots\) from `(t, q0x, q0y, p0x, p0y, m, alpha)`; fitness = error at horizon times.

Each scenario supplies its own `m` and `alpha` when evaluating fitness, so evolved laws can generalize across mass and coupling strength.

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
