# evophy

Small Clojure experiment: evolve symbolic expressions that behave like **conserved quantities** along a synthetic **harmonic-oscillator** trajectory. Ground-truth dynamics come from **[Emmy](https://github.com/mentat-collective/emmy/)** (Hamiltonian integration); a simple genetic program searches over trees built from `+`, `-`, `*`, `e/square`, coordinates `q` / `p`, and numeric constants.

Repository: [https://github.com/alishsan/evophy](https://github.com/alishsan/evophy)

## Requirements

- [Leiningen](https://leiningen.org/) 2.x
- JDK 11+ (tested with typical OpenJDK installs)

## Quick start

Clone the repo, then from the project root:

```bash
lein test
lein run
```

`lein run` integrates the reference oscillator, runs the GA for several generations, then prints the best individual’s **expression**, **fitness**, and a **result** map (`:mean` / `:variance` of the law over the trajectory samples).

## Tests

```bash
lein test
```

Covers trajectory generation and fitness behavior (e.g. flat literals vs state-dependent laws).

## Uberjar (optional)

```bash
lein uberjar
java -jar target/uberjar/evophy-*-standalone.jar
```

Version in the filename follows Leiningen’s default (`0.1.0-SNAPSHOT` until you release).

## License

This project is licensed under the **Eclipse Public License 2.0** (EPL-2.0). The full license text is in [`LICENSE`](LICENSE).
