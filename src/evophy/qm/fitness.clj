(ns evophy.qm.fitness
  "Fitness layer for the Schrodinger branch.

  Scores what a search proposes -- trial energies, trial wavefunctions, or
  whole *strategies* -- against the ground-truth oracle in
  `evophy.qm.schrodinger`. This is the QM analogue of the classical branch's
  fitness code: it keeps the same 'approximately right, and cheap counts too'
  spirit rather than exact-or-nothing. Concretely that means

    - partial credit: relative error is mapped through a smooth, bounded
      reward `1/(1 + (rel/tol)^2)` (the same shape the classical branch uses),
      so a trial that is close still scores well and search has a gradient to
      climb;
    - an oracle-free physics objective: the Rayleigh quotient
      `<psi|H|psi>/<psi|psi>` is a rigorous *upper bound* on the ground-state
      energy, so a trial wavefunction can be scored on physics alone (the
      oracle only validates), and excited states are reached by projecting out
      lower ones;
    - accuracy-vs-cost as two objectives, with a Pareto helper so search can
      trade a cheaper strategy against a more accurate one instead of
      collapsing everything to a single number too early.

  Everything here works on the same uniform grid the oracle uses, so trial
  wavefunctions and oracle wavefunctions are directly comparable."
  (:require [evophy.qm.schrodinger :as s]))

;; ---------------------------------------------------------------------------
;; Grid inner products (trapezoidal, matching the oracle's normalization)
;; ---------------------------------------------------------------------------

(defn sample
  "Sample a wavefunction fn `f` onto the oracle grid `xs`, returning a vector."
  [f xs]
  (mapv (fn [x] (double (f x))) xs))

(defn- trap-weight
  "Trapezoidal quadrature weight for grid index i of n points (endpoints 1/2)."
  ^double [^long i ^long n]
  (if (or (zero? i) (= i (dec n))) 0.5 1.0))

(defn inner-product
  "Trapezoidal <psi1|psi2> = integral psi1(x) psi2(x) dx over the grid.
  psi1/psi2 are vectors sampled on the grid described by `params`."
  ^double [params psi1 psi2]
  (let [h (s/step-size params)
        n (count psi1)]
    (loop [i 0 s 0.0]
      (if (< i n)
        (recur (inc i)
               (+ s (* (trap-weight i n)
                       (double (nth psi1 i)) (double (nth psi2 i)) h)))
        s))))

(defn norm2
  "<psi|psi> on the grid."
  ^double [params psi]
  (inner-product params psi psi))

(defn normalize
  "Scale psi to unit L2 norm on the grid. Returns a vector; a zero vector is
  returned unchanged (nothing to normalize)."
  [params psi]
  (let [nrm (Math/sqrt (max (norm2 params psi) 0.0))]
    (if (< nrm 1e-300)
      (vec psi)
      (mapv (fn [v] (/ (double v) nrm)) psi))))

;; ---------------------------------------------------------------------------
;; Energy scoring (partial credit)
;; ---------------------------------------------------------------------------

(def ^:const default-energy-tol
  "Relative-error scale at which an energy score falls to 1/2. Trials an order
  of magnitude closer score ~1; an order of magnitude worse score ~0.01."
  1e-2)

(defn energy-score
  "Partial-credit score in (0,1] for a trial energy against a reference.
  Uses `1/(1 + (rel/tol)^2)` on the relative error, so score = 1 at exact,
  1/2 at rel = tol, and decays smoothly (never hits 0, keeping a gradient).
  The reference is compared by magnitude, robust to E<0 bound states."
  (^double [^double trial-E ^double exact-E]
   (energy-score trial-E exact-E default-energy-tol))
  (^double [^double trial-E ^double exact-E ^double tol]
   (let [denom (max (Math/abs exact-E) 1e-12)
         rel   (/ (Math/abs (- trial-E exact-E)) denom)
         r     (/ rel tol)]
     (/ 1.0 (+ 1.0 (* r r))))))

(defn spectrum-score
  "Score a trial spectrum (vector of energies) against a reference spectrum.
  Trials and references are each sorted ascending and matched level-by-level;
  every reference level that has no trial counterpart scores 0, and surplus
  trial energies beyond the reference count are penalized the same way. This
  rewards getting both the values *and* the count right.

  Returns {:score mean-in-[0,1], :per-level [...], :n-ref, :n-trial}."
  [trial-energies exact-energies & {:keys [tol] :or {tol default-energy-tol}}]
  (let [trial (vec (sort (map double trial-energies)))
        exact (vec (sort (map double exact-energies)))
        nref  (count exact)
        ntri  (count trial)
        ;; one slot per reference level, plus one 0-slot per surplus trial
        per   (concat
                (map-indexed
                  (fn [i e]
                    (if (< i ntri) (energy-score (nth trial i) e tol) 0.0))
                  exact)
                (repeat (max 0 (- ntri nref)) 0.0))
        per   (vec per)]
    {:score (if (seq per) (/ (reduce + per) (count per)) 0.0)
     :per-level per
     :n-ref nref
     :n-trial ntri}))

;; ---------------------------------------------------------------------------
;; Wavefunction scoring (fidelity / overlap)
;; ---------------------------------------------------------------------------

(defn overlap
  "Fidelity |<psi-trial|psi-ref>| / (||psi-trial|| ||psi-ref||) in [0,1].
  Sign/phase-independent, so a wavefunction and its negation both score 1.
  1 = identical shape, 0 = orthogonal."
  ^double [params psi-trial psi-ref]
  (let [nt (Math/sqrt (max (norm2 params psi-trial) 0.0))
        nr (Math/sqrt (max (norm2 params psi-ref) 0.0))]
    (if (or (< nt 1e-300) (< nr 1e-300))
      0.0
      (/ (Math/abs (inner-product params psi-trial psi-ref)) (* nt nr)))))

(defn wavefn-score
  "Score a trial wavefunction against an oracle reference wavefunction.
  `psi-trial` may be a fn of x or a vector already sampled on the grid;
  `psi-ref` is a grid vector (e.g. the oracle's :psi). Returns the overlap
  fidelity in [0,1]."
  ^double [params psi-trial psi-ref]
  (let [xs (s/grid params)
        pt (if (fn? psi-trial) (sample psi-trial xs) psi-trial)]
    (overlap params pt psi-ref)))

;; ---------------------------------------------------------------------------
;; Variational objective (oracle-free physics)
;; ---------------------------------------------------------------------------

(defn rayleigh-quotient
  "Variational energy E[psi] = <psi|H|psi> / <psi|psi> for
  H = -(hbar^2/2m) d^2/dx^2 + V(x), evaluated on the grid.

  Uses the gradient form of the kinetic term,
  <psi|T|psi> = (hbar^2/2m) integral (psi')^2 dx, discretized with
  forward differences (equivalent to psi^T K psi for the standard tridiagonal
  Laplacian and always >= 0). Dirichlet boundaries (psi -> 0 just outside the
  grid) are assumed, as for any bound state, and are enforced explicitly via
  two boundary terms (psi[0]-0)^2/h and (psi[n-1]-0)^2/h -- without them a
  trial that fails to decay at the grid edge (e.g. a nonzero constant on the
  infinite-well benchmark, whose domain IS the hard wall) would dodge the
  kinetic cost of that discontinuity and can spuriously undercut the true
  ground energy. For potentials with padding beyond where a good trial decays
  (harmonic, finite/Woods-Saxon wells) these two terms are negligible. The
  result is a rigorous upper bound on the ground-state energy of V; the
  tighter (lower) it is, the better the trial.

  `psi` is a fn of x or a grid vector."
  ^double [params psi]
  (let [{:keys [hbar m V]} params
        xs  (s/grid params)
        p   (if (fn? psi) (sample psi xs) psi)
        h   (s/step-size params)
        n   (count p)
        kin-c (/ (* hbar hbar) (* 2.0 m))
        ;; kinetic: (hbar^2/2m) * sum ((p[i+1]-p[i])^2 / h), plus the two
        ;; Dirichlet boundary terms treating psi as 0 just outside the grid.
        kin (loop [i 0 s 0.0]
              (if (< i (dec n))
                (let [d (- (double (nth p (inc i))) (double (nth p i)))]
                  (recur (inc i) (+ s (/ (* d d) h))))
                (let [p0 (double (nth p 0))
                      pn (double (nth p (dec n)))]
                  (* kin-c (+ s (/ (* p0 p0) h) (/ (* pn pn) h))))))
        ;; potential + norm via trapezoidal
        [pot nrm] (loop [i 0 pv 0.0 nv 0.0]
                    (if (< i n)
                      (let [w  (trap-weight i n)
                            pi (double (nth p i))
                            vi (double (V (nth xs i)))]
                        (recur (inc i)
                               (+ pv (* w vi pi pi h))
                               (+ nv (* w pi pi h))))
                      [pv nv]))]
    (if (< nrm 1e-300)
      Double/POSITIVE_INFINITY
      (/ (+ kin pot) nrm))))

(defn project-out
  "Gram-Schmidt: remove the components of `psi` along each reference state in
  `refs` (grid vectors). Returns the residual grid vector, orthogonal to every
  ref. Used to turn the ground-state variational bound into a bound on an
  excited level: after projecting out the k lowest states, E[psi_perp] bounds
  the (k+1)-th eigenvalue."
  [params psi refs]
  (let [xs (s/grid params)
        p  (if (fn? psi) (sample psi xs) (vec psi))]
    (reduce
      (fn [cur ref]
        (let [rn (norm2 params ref)]
          (if (< rn 1e-300)
            cur
            (let [c (/ (inner-product params cur ref) rn)]
              (mapv (fn [a b] (- (double a) (* c (double b)))) cur ref)))))
      p refs)))

(defn variational-fitness
  "Oracle-free fitness for a trial ground-state wavefunction: the negative
  Rayleigh quotient, so that *higher is better* like the other scores here
  (search maximizes). Lower variational energy -> higher fitness.

  With `:below refs`, project the trial orthogonal to the given lower states
  first, yielding a bound on the next level up."
  [params psi & {:keys [below]}]
  (let [p (if (seq below) (project-out params psi below) psi)]
    {:energy (rayleigh-quotient params p)
     :fitness (- (rayleigh-quotient params p))}))

;; ---------------------------------------------------------------------------
;; Cost + Pareto (accuracy vs. work)
;; ---------------------------------------------------------------------------

(defn expr-cost
  "Structural cost of an expression tree: node count. A cheap proxy for the
  'work' / complexity of a trial-psi expression, so search can prefer simple
  forms (the classical branch's complexity penalty in the same spirit)."
  ^long [expr]
  (if (sequential? expr)
    (reduce + 1 (map expr-cost expr))
    1))

(defn scalar-fitness
  "Collapse accuracy and cost into one number when a search wants a scalar:
  accuracy discounted by cost, `accuracy / (1 + cost-weight * cost)`.
  Keep `cost-weight` small; use the Pareto helpers below to avoid collapsing
  prematurely."
  [{:keys [accuracy cost]} & {:keys [cost-weight] :or {cost-weight 0.0}}]
  (/ (double accuracy) (+ 1.0 (* (double cost-weight) (double (or cost 0))))))

(defn pareto-dominates?
  "True if candidate a Pareto-dominates b: a is no worse on accuracy (higher)
  and no worse on cost (lower), and strictly better on at least one axis.
  Each candidate is a map with :accuracy and :cost."
  [a b]
  (let [aa (double (:accuracy a)) ab (double (:accuracy b))
        ca (double (:cost a))     cb (double (:cost b))]
    (and (>= aa ab) (<= ca cb)
         (or (> aa ab) (< ca cb)))))

(defn pareto-front
  "Non-dominated subset of `candidates` (each a map with :accuracy and :cost).
  Preserves input order among the survivors."
  [candidates]
  (let [v (vec candidates)]
    (filterv
      (fn [c] (not-any? (fn [o] (and (not (identical? o c))
                                     (pareto-dominates? o c)))
                        v))
      v)))

;; ---------------------------------------------------------------------------
;; Strategy scoring (run a procedure, score its output against the oracle)
;; ---------------------------------------------------------------------------

(defn oracle-spectrum
  "Reference bound-state energies for a benchmark map, straight from the oracle.
  Convenience wrapper over `evophy.qm.schrodinger/spectrum`."
  [bench & opts]
  (mapv :energy (apply s/spectrum bench opts)))

(defn score-strategy
  "Score a *strategy*: a fn that, given the params, returns a proposed spectrum.
  The strategy may return either a bare vector of energies, or a map
  {:energies [...] :cost <n>} to report its own work (integration count, etc.).

  Returns {:accuracy (spectrum score in [0,1]), :cost, :fitness (scalar),
  :detail (per-level breakdown)}. `ref-energies` is the oracle truth (get it
  from `oracle-spectrum`)."
  [params strategy ref-energies
   & {:keys [tol cost-weight] :or {tol default-energy-tol cost-weight 0.0}}]
  (let [out       (strategy params)
        energies  (if (map? out) (:energies out) out)
        cost      (if (map? out) (or (:cost out) 0) 0)
        detail    (spectrum-score energies ref-energies :tol tol)
        accuracy  (:score detail)]
    {:accuracy accuracy
     :cost cost
     :fitness (scalar-fitness {:accuracy accuracy :cost cost} :cost-weight cost-weight)
     :detail detail}))

;; ---------------------------------------------------------------------------
;; Demo
;; ---------------------------------------------------------------------------

(defn -main
  "Demonstrate the three scoring modes on the harmonic oscillator, whose
  answers we know: variational energy from a trial wavefunction, wavefunction
  fidelity against the oracle, and end-to-end strategy scoring."
  [& _]
  (let [bench  (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
        params (:params bench)
        psi0   (:psi (first (s/solve params :n-levels 1 :E-hi 20.0)))]
    (println "Harmonic oscillator (omega=1, hbar=m=1): E_0 = 0.5, E_n = n + 1/2")
    (println)
    (println "1) Variational energy E[psi] = <psi|H|psi>/<psi|psi> (upper bound on E_0):")
    (doseq [[label a] [["exact  exp(-x^2/2)" 0.5]
                       ["wide   exp(-x^2/4)" 0.25]
                       ["narrow exp(-x^2)"   1.0]]]
      (let [e (rayleigh-quotient params (fn [x] (Math/exp (* (- a) x x))))]
        (println (format "   %-20s E = %.5f   fitness = %.5f" label e (- e)))))
    (println)
    (println "2) Wavefunction fidelity |<trial|oracle_0>| (1 = identical shape):")
    (doseq [[label f] [["exact  exp(-x^2/2)" (fn [x] (Math/exp (* -0.5 x x)))]
                       ["wrong  exp(-x^2/8)" (fn [x] (Math/exp (* -0.125 x x)))]]]
      (println (format "   %-20s overlap = %.5f" label (wavefn-score params f psi0))))
    (println)
    (println "3) Strategy scoring against the oracle spectrum (accuracy in [0,1]):")
    (let [ref (oracle-spectrum bench :n-levels 4 :E-hi 20.0)]
      (doseq [[label strat] [["oracle truth" (fn [_] ref)]
                             ["+0.3 shift"   (fn [_] (mapv #(+ % 0.3) ref))]
                             ["only 2 levels" (fn [_] (vec (take 2 ref)))]]]
        (let [{:keys [accuracy]} (score-strategy params strat ref)]
          (println (format "   %-16s accuracy = %.5f" label accuracy)))))))
