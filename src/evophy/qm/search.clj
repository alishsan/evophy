(ns evophy.qm.search
  "GP search layer for the Schrodinger branch (roadmap item 3).

  Evolves trial-psi expression trees over `x` and scores them with the
  oracle-*free* variational objective from `evophy.qm.fitness`
  (`rayleigh-quotient`): E[psi] = <psi|H|psi>/<psi|psi> is a rigorous upper
  bound on the ground-state energy for ANY normalizable trial, so search
  needs no reference answer to climb toward truth -- only physics. Excited
  states are reached the same way the fitness layer intends: sequential
  Gram-Schmidt deflation (`project-out`) against the psis already found,
  turning each new variational bound into a bound on the next level up.

  This is a fresh pipeline, not a reuse of `evophy.core`'s Kepler machinery
  (see the qm README): a single-slot genome `{:expr <s-expr over x>}`, a
  small generic arithmetic+transcendental primitive set (no orbital
  primitives, no phase-space charts), and a plain truncation-elitism GP loop
  in the same spirit as the classical branch's `evolve-generation` (sort by
  fitness, keep elites, breed the rest by mutation/crossover, top up with
  random immigrants) -- reimplemented from scratch against this domain's own
  representation."
  (:require [evophy.qm.schrodinger :as s]
            [evophy.qm.fitness :as f]))

;; ---------------------------------------------------------------------------
;; Primitive set: generic real-valued ops over a single variable `x`, plus the
;; potential V(x) itself (exactly the known input to the equation -- search
;; still has to discover how psi relates to it, not be handed the answer).
;; Every op is domain-total (never NaN/Infinite) so mutation/crossover never
;; wander into dead subtrees that silently poison the whole expression.
;; ---------------------------------------------------------------------------

(defn- qexp ^double [a] (Math/exp (max -40.0 (min 40.0 (double a)))))
(defn- qneg ^double [a] (- (double a)))
(defn- qsquare ^double [a] (let [x (double a)] (* x x)))
(defn- qabs ^double [a] (Math/abs (double a)))
(defn- qsqrt ^double [a] (Math/sqrt (Math/abs (double a))))
(defn- qsin ^double [a] (Math/sin (double a)))
(defn- qcos ^double [a] (Math/cos (double a)))
(defn- qdiv ^double [a b]
  (let [x (double a) y (double b)]
    (/ x (if (< (Math/abs y) 1e-9) (if (neg? y) -1e-9 1e-9) y))))
(defn- qsafe ^double [a]
  (let [d (double a)] (if (Double/isFinite d) d 0.0)))

(def ^:private unary-ops '#{qexp qneg qsquare qabs qsqrt qsin qcos V})
(def ^:private binary-ops '[+ - * qdiv])
(def ^:private all-ops (into (vec unary-ops) binary-ops))

(def ^:private default-max-depth 5)
(def ^:private max-expr-size 45)

(defn- random-constant [] (- (* 4.0 (rand)) 2.0))

(defn- random-terminal []
  (if (< (rand) 0.5) 'x (random-constant)))

(defn random-expr
  "A random s-expression over `x` (and the potential `V`), depth-bounded."
  [depth]
  (if (or (<= depth 0) (< (rand) 0.3))
    (random-terminal)
    (let [op (rand-nth all-ops)]
      (if (unary-ops op)
        (list op (random-expr (dec depth)))
        (list op (random-expr (dec depth)) (random-expr (dec depth)))))))

;; ---------------------------------------------------------------------------
;; Tree editing: index subexpressions in a fixed pre-order (argument
;; positions only -- operator symbols are never selected, so mutation and
;; crossover can't corrupt arity) for mutation/crossover.
;; ---------------------------------------------------------------------------

(defn- count-subexprs ^long [expr]
  (if (sequential? expr)
    (inc (reduce + 0 (map count-subexprs (rest expr))))
    1))

(defn- subexpr-at [expr idx]
  (if (zero? idx)
    expr
    (let [[_ & args] expr]
      (loop [i (dec (long idx)) remaining args]
        (let [a (first remaining) sz (count-subexprs a)]
          (if (< i sz)
            (subexpr-at a i)
            (recur (- i sz) (rest remaining))))))))

(defn- replace-subexpr-at [expr idx replacement]
  (if (zero? idx)
    replacement
    (let [[op & args] expr]
      (loop [i (dec (long idx)) pre [] remaining args]
        (let [a (first remaining) sz (count-subexprs a)]
          (if (< i sz)
            (cons op (concat pre [(replace-subexpr-at a i replacement)] (rest remaining)))
            (recur (- i sz) (conj pre a) (rest remaining))))))))

(defn- bounded
  "Fall back to a fresh small random tree if editing produced a bloated one."
  [expr]
  (if (<= (count-subexprs expr) max-expr-size)
    expr
    (random-expr 2)))

(defn- jitter [c]
  (+ (double c) (* 0.5 (- (rand) 0.5) (max 1.0 (Math/abs (double c))))))

(defn mutate-expr
  "Pick one random subtree and replace it: jitter it if it's a constant,
  otherwise swap in a fresh random subtree."
  [expr]
  (let [n (count-subexprs expr)
        idx (rand-int n)
        target (subexpr-at expr idx)
        replacement (if (and (number? target) (< (rand) 0.5))
                      (jitter target)
                      (random-expr (max 1 (rand-int default-max-depth))))]
    (bounded (replace-subexpr-at expr idx replacement))))

(defn crossover-expr
  "Subtree-swap crossover: graft a random subtree of `b` into a random spot
  in `a`."
  [a b]
  (bounded
    (replace-subexpr-at a (rand-int (count-subexprs a))
                         (subexpr-at b (rand-int (count-subexprs b))))))

;; ---------------------------------------------------------------------------
;; Compilation: expr -> (fn [x V] double), memoized by printed form. `V` is
;; passed in as an argument (not closed over) so the same compiled fn is
;; valid for any potential -- exactly the trick the classical branch's
;; `compile-state-fn` uses to stay generic across initial conditions.
;; ---------------------------------------------------------------------------

(def ^:private ^java.util.concurrent.ConcurrentHashMap compile-cache
  (java.util.concurrent.ConcurrentHashMap.))

(defn compile-expr
  "Compile a trial-psi expr to a fn [x V] -> double, where V is itself a
  fn of x (the potential). Compiled fns are cached by printed form."
  [expr]
  (let [key (pr-str expr)]
    (or (.get compile-cache key)
        (let [f (binding [*ns* (the-ns 'evophy.qm.search)]
                  (eval `(fn [~'x ~'V] (qsafe ~expr))))]
          (.put compile-cache key f)
          f))))

;; ---------------------------------------------------------------------------
;; Individuals + fitness
;; ---------------------------------------------------------------------------

(def ^:private seed-exprs
  "Textbook variational ansatz families, seeded into every search so it never
  has to reinvent decay from scratch: Gaussian, odd-Gaussian, exponential
  (Yukawa-like), Lorentzian."
  '[(qexp (qneg (qsquare x)))
    (* x (qexp (qneg (qsquare x))))
    (qexp (qneg (qabs x)))
    (qdiv 1.0 (+ 1.0 (qsquare x)))])

(defn- seed-individuals [] (mapv (fn [e] {:expr e}) seed-exprs))

(defn random-individual [] {:expr (random-expr default-max-depth)})

(defn evaluate-individual
  "Score `ind` against `params`'s Hamiltonian via the oracle-free variational
  objective, optionally projected orthogonal to `below` (grid vectors of
  lower states already found). `:fitness` is the Rayleigh-quotient fitness
  discounted by a small parsimony penalty on expression size."
  [params ind & {:keys [below cost-weight] :or {below [] cost-weight 5.0e-4}}]
  (let [fc (compile-expr (:expr ind))
        Vf (:V params)
        psi (fn [x] (fc x Vf))
        {:keys [energy fitness]} (f/variational-fitness params psi :below below)
        cost (f/expr-cost (:expr ind))
        adj (if (Double/isFinite fitness)
              (- fitness (* (double cost-weight) cost))
              Double/NEGATIVE_INFINITY)]
    (assoc ind :energy energy :fitness adj :cost cost)))

;; ---------------------------------------------------------------------------
;; Evolutionary loop: truncation elitism + mutation/crossover + immigrants,
;; the same shape as the classical branch's `evolve-generation` but
;; reimplemented against this domain's single-expr genome.
;; ---------------------------------------------------------------------------

(defn- breed [elites]
  (let [p1 (rand-nth elites)]
    {:expr (if (< (rand) 0.65)
             (crossover-expr (:expr p1) (:expr (rand-nth elites)))
             (mutate-expr (:expr p1)))}))

(defn evolve-level
  "GP search for one bound state of `params`'s Hamiltonian. With `:below`
  (grid vectors of already-found lower states), the trial is projected
  orthogonal to them first, so this bounds the next level up instead of the
  ground state. Returns the best individual seen: {:expr :energy :fitness
  :cost}."
  [params & {:keys [below population-size generations elite-frac immigrant-frac
                     cost-weight]
             :or {below [] population-size 150 generations 60
                  elite-frac 0.15 immigrant-frac 0.12 cost-weight 5.0e-4}}]
  (let [elite-n (max 2 (long (* (double elite-frac) population-size)))
        imm-n (max 1 (long (* (double immigrant-frac) population-size)))
        seeds (seed-individuals)
        init (into seeds (repeatedly (max 0 (- population-size (count seeds)))
                                      random-individual))]
    (loop [gen 0 population init best nil]
      (if (>= gen generations)
        best
        (let [scored (->> population
                           (mapv #(evaluate-individual params % :below below
                                                        :cost-weight cost-weight))
                           (sort-by :fitness >))
              gen-best (first scored)
              best' (if (or (nil? best) (> (double (:fitness gen-best)) (double (:fitness best))))
                      gen-best best)
              top (vec (take elite-n scored))
              elite-pool (if (some #(= (:expr %) (:expr best')) top) top (into [best'] top))
              n-children (max 0 (- population-size imm-n (count elite-pool)))
              children (vec (repeatedly n-children #(breed elite-pool)))
              immigrants (vec (repeatedly imm-n random-individual))]
          (recur (inc gen) (into elite-pool (concat children immigrants)) best'))))))

(defn evolve-spectrum
  "Evolve the first `n-levels` bound states in ascending order via sequential
  deflation: each level's search runs orthogonal to the (grid, normalized)
  psis of every level already found. Returns a vector of {:level :expr
  :energy :fitness :cost :psi}, `opts` forwarded to `evolve-level`."
  [params n-levels & opts]
  (let [xs (s/grid params)
        Vf (:V params)]
    (loop [lvl 0 below [] results []]
      (if (>= lvl (long n-levels))
        results
        (let [best (apply evolve-level params :below below opts)
              fc (compile-expr (:expr best))
              psi-raw (f/sample (fn [x] (fc x Vf)) xs)
              psi-proj (if (seq below) (f/project-out params psi-raw below) psi-raw)
              psi-norm (f/normalize params psi-proj)]
          (recur (inc lvl) (conj below psi-norm)
                 (conj results (assoc (select-keys best [:expr :energy :fitness :cost])
                                       :level lvl :psi psi-norm))))))))

;; ---------------------------------------------------------------------------
;; Demo
;; ---------------------------------------------------------------------------

(defn -main
  "Search the harmonic oscillator for its first three levels from scratch --
  no closed-form answer given to the search, only V(x) and the variational
  principle -- and compare the discovered energies to the known spectrum."
  [& _]
  (let [bench (s/harmonic {:omega 1.0 :x-max 10.0 :n 1001})
        params (:params bench)
        exact (:exact bench)]
    (println "Harmonic oscillator (omega=1, hbar=m=1): searching for E_0, E_1, E_2 (exact 0.5, 1.5, 2.5)")
    (println "No spectrum given to the search -- only V(x) and the variational principle.")
    (println)
    (doseq [{:keys [level expr energy]} (evolve-spectrum params 3)]
      (println (format "  level %d: E = %.5f  (exact %.5f)  psi = %s"
                        level energy (exact level) (pr-str expr))))))
