(ns evophy.qm.search-test
  "Validate the QM search layer. Tree-editing helpers get exact structural
  checks; the GP loop itself is checked against two properties that hold
  regardless of RNG draw, so the suite stays non-flaky:
    - the Rayleigh quotient is a rigorous upper bound on the true energy
      (never violated, however bad the trial), and
    - truncation elitism guarantees the returned best is never worse than
      the best seed present in generation 0 (with cost-weight 0, so ranking
      is exactly by energy)."
  (:require [clojure.test :refer [deftest is]]
            [evophy.qm.schrodinger :as s]
            [evophy.qm.search :as search]))

;; ---------------------------------------------------------------------------
;; Tree editing: exact, deterministic
;; ---------------------------------------------------------------------------

(deftest subexpr-indexing-round-trips
  (let [e '(+ x 1)]
    (is (= 3 (#'search/count-subexprs e)))
    (is (= e (#'search/subexpr-at e 0)))
    (is (= 'x (#'search/subexpr-at e 1)))
    (is (= 1 (#'search/subexpr-at e 2)))
    (is (= '(+ y 1) (#'search/replace-subexpr-at e 1 'y)))
    (is (= '(+ x 99) (#'search/replace-subexpr-at e 2 99)))
    (is (= 'z (#'search/replace-subexpr-at e 0 'z)))))

(deftest subexpr-indexing-nested
  (let [e '(* (+ x 1) x)]
    (is (= 5 (#'search/count-subexprs e)))
    (is (= '(+ x 1) (#'search/subexpr-at e 1)))
    (is (= 'x (#'search/subexpr-at e 2)))
    (is (= 1 (#'search/subexpr-at e 3)))
    (is (= 'x (#'search/subexpr-at e 4)))
    (is (= '(* (+ x 7) x) (#'search/replace-subexpr-at e 3 7)))
    (is (= '(* (+ x 1) 7) (#'search/replace-subexpr-at e 4 7)))))

(deftest mutate-and-crossover-stay-well-formed
  (dotimes [_ 25]
    (let [a (search/random-expr 4)
          b (search/random-expr 4)
          m (search/mutate-expr a)
          c (search/crossover-expr a b)]
      ;; a well-formed expr compiles and evaluates to a finite double
      (doseq [e [m c]]
        (let [fc (search/compile-expr e)]
          (is (Double/isFinite (fc 1.0 (fn [_] 0.0)))))))))

;; ---------------------------------------------------------------------------
;; Compilation + evaluation
;; ---------------------------------------------------------------------------

(deftest compile-expr-matches-hand-eval
  (let [fc (search/compile-expr '(qexp (qneg (qsquare x))))]
    (is (< (Math/abs (- (fc 0.0 (fn [_] 0.0)) 1.0)) 1e-9))
    (is (< (Math/abs (- (fc 1.0 (fn [_] 0.0)) (Math/exp -1.0))) 1e-9))))

(deftest compile-expr-uses-passed-in-potential
  (let [fc (search/compile-expr '(V x))]
    (is (== 5.0 (fc 2.0 (fn [x] (* 2.5 x)))))
    (is (== 10.0 (fc 4.0 (fn [x] (* 2.5 x)))))))

(deftest evaluate-individual-recovers-known-energy
  (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
        ind (search/evaluate-individual params {:expr '(qexp (qneg (* 0.5 (qsquare x))))}
                                         :cost-weight 0.0)]
    (is (< (Math/abs (- 0.5 (:energy ind))) 5e-3))
    (is (< (Math/abs (- (- (:energy ind)) (:fitness ind))) 1e-9)
        "with cost-weight 0, fitness is exactly -energy")))

;; ---------------------------------------------------------------------------
;; The GP loop: two RNG-independent guarantees
;; ---------------------------------------------------------------------------

(deftest evolve-level-never-beats-the-variational-bound
  ;; The theorem holds for the continuous Hamiltonian; the discretized
  ;; kinetic term used by rayleigh-quotient can undershoot it by a hair
  ;; (fitness_test.clj's own upper-bound test uses the same 0.001 slack).
  (let [{:keys [params exact]} (s/harmonic {:omega 1.0 :x-max 10.0 :n 1001})
        best (search/evolve-level params :population-size 60 :generations 15
                                   :cost-weight 0.0)]
    (is (>= (:energy best) (- (exact 0) 1e-2))
        "Rayleigh quotient can never fall meaningfully below the true ground energy")))

(deftest evolve-level-is-never-worse-than-its-best-seed
  (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 10.0 :n 1001})
        seed-exprs '[(qexp (qneg (qsquare x)))
                     (* x (qexp (qneg (qsquare x))))
                     (qexp (qneg (qabs x)))
                     (qdiv 1.0 (+ 1.0 (qsquare x)))]
        seed-energies (map (fn [e] (:energy (search/evaluate-individual
                                              params {:expr e} :cost-weight 0.0)))
                            seed-exprs)
        best (search/evolve-level params :population-size 60 :generations 15
                                   :cost-weight 0.0)]
    (is (<= (:energy best) (+ 1e-9 (apply min seed-energies)))
        "elitism guarantees the final best is never worse than the best seed")))

(deftest evolve-spectrum-produces-ascending-levels
  (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 10.0 :n 601})
        results (search/evolve-spectrum params 2 :population-size 40 :generations 10
                                         :cost-weight 0.0)]
    (is (= [0 1] (mapv :level results)))
    (is (every? #(Double/isFinite (:energy %)) results))
    (is (= 601 (count (:psi (first results)))))))

(deftest evolve-spectrum-second-level-is-in-the-right-ballpark
  ;; Deflation projects against a *found* (not exact) ground state, so the
  ;; E_1 upper-bound theorem only holds approximately here -- generous
  ;; slack, just enough to catch a broken projection (e.g. level 1 landing
  ;; back on E_0) without being sensitive to search variance.
  (let [{:keys [params exact]} (s/harmonic {:omega 1.0 :x-max 10.0 :n 601})
        [_ lvl1] (search/evolve-spectrum params 2 :population-size 60 :generations 15
                                          :cost-weight 0.0)]
    (is (< (Math/abs (- (:energy lvl1) (exact 1))) 0.9))))
