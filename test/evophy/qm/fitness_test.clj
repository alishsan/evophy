(ns evophy.qm.fitness-test
  "Validate the fitness layer against the oracle. Two kinds of checks:
  algebraic properties of the scoring functions (partial credit shape,
  fidelity, Pareto dominance) and physics properties tied to the oracle
  (the Rayleigh quotient is a real variational upper bound; oracle
  eigenstates are orthonormal and self-score 1)."
  (:require [clojure.test :refer [deftest is testing]]
            [evophy.qm.schrodinger :as s]
            [evophy.qm.fitness :as f]))

;; ---------------------------------------------------------------------------
;; Energy scoring: partial credit
;; ---------------------------------------------------------------------------

(deftest energy-score-shape
  (testing "exact -> 1, tol away -> 1/2, monotone decreasing in error"
    (is (== 1.0 (f/energy-score 5.0 5.0)) "exact scores 1")
    (is (< (Math/abs (- 0.5 (f/energy-score (* 5.0 1.01) 5.0 1e-2))) 1e-9)
        "rel = tol scores 1/2")
    (is (> (f/energy-score 5.05 5.0) (f/energy-score 5.20 5.0))
        "closer trial scores higher")
    (is (< 0.0 (f/energy-score 50.0 5.0)) "always positive (keeps a gradient)"))
  (testing "works for negative (bound-state) energies via magnitude"
    (is (== 1.0 (f/energy-score -12.5 -12.5)))
    (is (> (f/energy-score -12.4 -12.5) (f/energy-score -11.0 -12.5)))))

(deftest spectrum-score-rewards-values-and-count
  (testing "perfect spectrum scores ~1"
    (let [d (f/spectrum-score [1.0 2.0 3.0] [1.0 2.0 3.0])]
      (is (< (Math/abs (- 1.0 (:score d))) 1e-9))))
  (testing "a missing level drags the mean down"
    (let [d (f/spectrum-score [1.0 2.0] [1.0 2.0 3.0])]
      (is (= 3 (:n-ref d)))
      (is (= 2 (:n-trial d)))
      (is (< (:score d) 0.7) "one of three levels unmatched -> <= 2/3")))
  (testing "surplus trial energies are penalized"
    (let [d (f/spectrum-score [1.0 2.0 3.0 4.0] [1.0 2.0 3.0])]
      (is (< (:score d) 1.0) "a spurious extra level lowers the score"))))

;; ---------------------------------------------------------------------------
;; Wavefunction fidelity against the oracle
;; ---------------------------------------------------------------------------

(deftest oracle-eigenstates-are-orthonormal
  (testing "self-overlap 1, cross-overlap ~0 for oracle harmonic states"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
          levels (s/solve params :n-levels 3 :E-hi 20.0)
          [g e2 _] (map :psi levels)]
      (is (< (Math/abs (- 1.0 (f/overlap params g g))) 1e-6) "ground self-overlap 1")
      (is (< (f/overlap params g e2) 1e-2) "distinct levels nearly orthogonal"))))

(deftest wavefn-score-matches-analytic-ground-state
  (testing "analytic Gaussian ground state overlaps the oracle ground state"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
          psi0 (:psi (first (s/solve params :n-levels 1 :E-hi 20.0)))
          gauss (fn [x] (Math/exp (* -0.5 x x)))]
      (is (> (f/wavefn-score params gauss psi0) 0.999)
          "exp(-x^2/2) is the harmonic ground state up to normalization"))))

;; ---------------------------------------------------------------------------
;; Variational objective (oracle-free physics)
;; ---------------------------------------------------------------------------

(deftest rayleigh-recovers-ground-energy
  (testing "E[exact ground state] ~ 0.5 for the harmonic oscillator"
    (let [params (:params (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001}))
          gauss  (fn [x] (Math/exp (* -0.5 x x)))]
      (is (< (Math/abs (- 0.5 (f/rayleigh-quotient params gauss))) 5e-3)))))

(deftest rayleigh-is-an-upper-bound
  (testing "any non-eigenstate trial gives E >= E_0"
    (let [params (:params (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001}))
          exact  (fn [x] (Math/exp (* -0.5 x x)))
          ;; wrong Gaussian width: variational energy must be strictly higher
          wide   (fn [x] (Math/exp (* -0.25 x x)))
          narrow (fn [x] (Math/exp (* -1.0 x x)))]
      (is (> (f/rayleigh-quotient params wide) (f/rayleigh-quotient params exact)))
      (is (> (f/rayleigh-quotient params narrow) (f/rayleigh-quotient params exact)))
      (is (>= (f/rayleigh-quotient params wide) 0.499)
          "upper bound never dips below the true ground energy"))))

(deftest projection-bounds-excited-level
  (testing "projecting out the ground state bounds the first excited energy"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001})
          levels (s/solve params :n-levels 2 :E-hi 20.0)
          psi0   (:psi (first levels))
          ;; a generic odd trial, not the exact first excited state
          trial  (fn [x] (* x (Math/exp (* -0.4 x x))))
          {:keys [energy]} (f/variational-fitness params trial :below [psi0])]
      (is (>= energy 1.49) "bound sits at/above the true E_1 = 1.5")
      (is (< energy 2.0)   "and a decent odd trial gets reasonably close"))))

(deftest variational-fitness-higher-is-better
  (testing "fitness increases as the trial energy decreases"
    (let [params (:params (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001}))
          good (f/variational-fitness params (fn [x] (Math/exp (* -0.5 x x))))
          bad  (f/variational-fitness params (fn [x] (Math/exp (* -0.2 x x))))]
      (is (> (:fitness good) (:fitness bad))))))

;; ---------------------------------------------------------------------------
;; Cost / Pareto
;; ---------------------------------------------------------------------------

(deftest expr-cost-counts-nodes
  (is (= 1 (f/expr-cost 'x)))
  (is (= 4 (f/expr-cost '(+ x 1))))          ;; +, x, 1 -> 3 nodes + 1 for the list
  (is (= 7 (f/expr-cost '(* (+ x 1) x)))))   ;; *, (+ x 1)=4, x -> 1+4+1 + 1 list

(deftest pareto-dominance
  (testing "dominance requires no-worse on both axes and better on one"
    (is (f/pareto-dominates? {:accuracy 0.9 :cost 2} {:accuracy 0.8 :cost 2}))
    (is (f/pareto-dominates? {:accuracy 0.9 :cost 1} {:accuracy 0.9 :cost 3}))
    (is (not (f/pareto-dominates? {:accuracy 0.9 :cost 3} {:accuracy 0.8 :cost 1}))
        "cheaper-but-worse does not dominate"))
  (testing "front keeps only non-dominated candidates"
    (let [cs [{:id :a :accuracy 0.9 :cost 3}
              {:id :b :accuracy 0.7 :cost 1}
              {:id :c :accuracy 0.6 :cost 2}]  ;; dominated by b
          front (set (map :id (f/pareto-front cs)))]
      (is (contains? front :a))
      (is (contains? front :b))
      (is (not (contains? front :c))))))

;; ---------------------------------------------------------------------------
;; Strategy scoring end-to-end
;; ---------------------------------------------------------------------------

(deftest score-strategy-against-oracle
  (testing "a strategy that returns the oracle spectrum scores ~1"
    (let [bench (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
          ref   (f/oracle-spectrum bench :n-levels 4 :E-hi 20.0)
          truth-strategy (fn [_] ref)
          off-strategy   (fn [_] (mapv #(+ % 0.3) ref))
          truth (f/score-strategy (:params bench) truth-strategy ref)
          off   (f/score-strategy (:params bench) off-strategy ref)]
      (is (> (:accuracy truth) 0.999) "exact reproduction scores ~1")
      (is (< (:accuracy off) (:accuracy truth)) "shifted energies score worse")))
  (testing "cost weighting discounts an expensive strategy"
    (let [bench (s/harmonic {:omega 1.0 :x-max 12.0 :n 2001})
          ref   (f/oracle-spectrum bench :n-levels 3 :E-hi 20.0)
          expensive (fn [_] {:energies ref :cost 100})
          scored (f/score-strategy (:params bench) expensive ref :cost-weight 0.01)]
      (is (< (:fitness scored) (:accuracy scored))
          "fitness is discounted below raw accuracy by the cost term"))))
